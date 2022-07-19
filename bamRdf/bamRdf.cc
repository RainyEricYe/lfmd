#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/GenomicRegionCollection.h"
using namespace SeqLib;
using namespace std;

#define PROGRAM "bamRdf"
#define VERSION "v1.4"
#define AUTHORS "yerui"
#define CONTACT "yerui@connect.hku.hk"
#define COMMENT "(cluster read pair by positon)"

class Fam;
typedef map<GenomicRegion, vector<BamRecord> > mRegBrV;
typedef map<GenomicRegion, vector<Fam> > mRegFamV;
typedef unsigned long ulong;

class Fam {
    public:
        ulong fh;  // filehandle id
        ulong read_n;
        string chr_start_end;

        Fam():fh(0), read_n(0), chr_start_end("") {};
        ~Fam(){};
};

bool add2mapRegionBamRecordV(const BamRecord &ra, const BamRecord &rb, mRegBrV &h, const GRC &grc, bool mitoMode);
inline int32_t true_start(const BamRecord &a);
inline int32_t true_end(const BamRecord &a);
inline string itoa(ulong &n) { ostringstream s; s << n; return s.str(); }

ulong sum_reads_n( vector<Fam> &vFam );
void show( SeqLib::BamWriter &w, ostream &outRdf, vector<SeqLib::BamReader> &vecR, ifstream *vecF, vector<Fam> &vFam, bool ok );
void add2mapRegionFamV( ifstream *vecF, ulong &i, mRegFamV &rfmap, SeqLib::BamHeader &head );
void mergeBlockFiles(SeqLib::BamWriter &w, ofstream &outRdf, string &outBam, ulong &part, BamHeader &head, ulong &rdfSizeCut );
void printBlockFile(string &outBam, ulong &part, mRegBrV &h, const BamHeader &head );

void usage()
{
    cout <<
        "Program: " << PROGRAM << " " << COMMENT << "\n"
        "Version: " << VERSION << "\n"
        "Authors: " << AUTHORS << "\n"
        "Contact: " << CONTACT << "\n\n"

        "Usage: " << PROGRAM << " -i in.bam\n\n"

        "       -t [s]    target region file (bed)\n"
        "       -i [s]    input original bam file which sorted by Qname\n"
        "       -o [s]    output bam file [out.bam]\n"
        "       -d [i]    discard reads family which size < [int]  [0]\n"
        "       -b [i]    memory control: N pairs of reads as a block [3000000]\n"
        "       -m        mitochondrial mode. allow read start < 0 [false]\n"
        "       -n [i]    limit number of read pairs [0, unlimited]\n"
        "       -v        version\n"
        "       -h        help\n"
        "\n";
}

int main( int argc, char** argv ) {

    string target_f(""), inBam(""), outBam("out.bam");
    ulong rdfSizeCut(0);
    ulong part(0);
    ulong ReadPairNum(0);
    ulong cntReadPairNum(0);

    ulong blockSize(3000000);
    ulong cnt(0); // pairs of read
    bool mitoMode(false);

    int c;
    while ( (c=getopt(argc,argv,"t:i:o:r:d:b:n:mvh")) != -1 ) {
        switch (c) {
            case 't': target_f = optarg;          break;
            case 'i': inBam    = optarg;          break;
            case 'o': outBam   = optarg;          break;
            case 'd': rdfSizeCut = atoi(optarg);  break;
            case 'b': blockSize = atoi(optarg);   break;
            case 'n': ReadPairNum = atoi(optarg); break;
            case 'm': mitoMode = true;            break;
            case 'v': cerr << VERSION << endl;    exit(1);
            case 'h':
            default:  usage();                    exit(1);
        }
    }

    if ( blockSize < 50000 )  blockSize = 50000;
    if ( inBam.size() == 0 ) cerr << "-i in.bam is needed" << endl, usage(), exit(1);

    SeqLib::BamReader r;
    SeqLib::BamWriter w(SeqLib::BAM);

    string outRdf_f = outBam + ".rdf";
    ofstream outRdf(outRdf_f.c_str());

    GRC grc;
    mRegBrV h;

    if ( !r.Open(inBam) )    cerr << "cannot open file: " << inBam    << endl, exit(1);
    if ( !w.Open(outBam) )   cerr << "cannot open file: " << outBam   << endl, exit(1);
    if ( !outRdf.is_open() ) cerr << "cannot open file: " << outRdf_f << endl, exit(1);

    SeqLib::BamHeader head = r.Header();
    w.SetHeader(head);
    w.WriteHeader();

    if ( target_f != "" ) {
        grc.ReadBED(target_f, head);
        grc.CreateTreeMap();
    }

    BamRecord ra, rb;
    r.GetNextRecord(ra);

    while ( r.GetNextRecord(rb) ) {
        if ( rb.SecondaryFlag() || rb.Interchromosomal() ) continue;

        if ( ra.Qname() == rb.Qname() ) {
            if ( add2mapRegionBamRecordV(ra, rb, h, grc, mitoMode) ) {
                cnt++;
                if ( cnt == blockSize ) {
                    printBlockFile(outBam, part, h, head );
                    h.clear();
                    part++;
                    cnt = 0;
                }

                cntReadPairNum++;
                if ( cntReadPairNum == ReadPairNum ) break;
            }

            if ( !r.GetNextRecord(ra) ) break;
        }
        else {
            ra = rb;
        }
    }
    r.Close();

    if ( !h.empty() ) {
        printBlockFile(outBam, part, h, head );
        h.clear();
        part++;
    }

    if ( part == 1 ) {
        cerr << "only has one block" << endl;
        w.Close();
        outRdf.close();

        string cmd = "mv -f " + outBam + ".0 " + outBam + " && mv -f " + outBam + ".0.rdf " + outBam + ".rdf";
        system( cmd.c_str() );
    }
    else {
        cerr << "merge " << part << " blocks..." << endl;
        mergeBlockFiles(w, outRdf, outBam, part, head, rdfSizeCut);
    }

    exit(0);
}

//~~~~~~~~~~~~~~

bool add2mapRegionBamRecordV(const BamRecord &ra, const BamRecord &rb, mRegBrV &h, const GRC &grc, bool mitoMode )
{
    if ( ra.ChrID() != rb.ChrID() or !ra.MappedFlag() or !rb.MappedFlag() )
        return false;

    GenomicRegion gra = ra.AsGenomicRegion();
    GenomicRegion grb = rb.AsGenomicRegion();

    if ( grc.size() > 0 && !grc.CountOverlaps(gra) && !grc.CountOverlaps(grb) )
        return false;

    // mitochondrial mode which allow true start <= 0
    if ( !mitoMode && (true_start(ra) <= 0 || true_start(rb) <= 0) )
        return false;

    // proper pair and orientation
    if ( ra.ReverseFlag() && !rb.ReverseFlag()
            && true_start(rb) <= true_start(ra)
       ) {

        GenomicRegion key( ra.ChrID(), true_start(rb), true_end(ra) );
        h[key].push_back(rb);
        h[key].push_back(ra);
    }
    else if ( rb.ReverseFlag() && !ra.ReverseFlag()
                 && true_start(ra) <= true_start(rb)
            ) {

        GenomicRegion key( ra.ChrID(), true_start(ra), true_end(rb) );
        h[key].push_back(ra);
        h[key].push_back(rb);
    }
    else {
        return false;
    }
}

inline int32_t true_start( const BamRecord &a )
{
    CigarField cf = a.GetCigar().front();

    if ( cf.Type() == 'S' )
        return a.Position() - cf.Length();
    else
        return a.Position();
}

inline int32_t true_end( const BamRecord &a )
{
    return true_start(a) + a.Length() - 1;
}

ulong sum_reads_n( vector<Fam> &vFam )
{
    ulong n(0);
    for ( auto& fam : vFam ) {
        n += fam.read_n;
    }

    return n;
}

void show( SeqLib::BamWriter &w, ostream &outRdf, vector<SeqLib::BamReader> &vecR, ifstream *vecF, vector<Fam> &vFam, bool ok )
{

    ulong total_read_n(0);

    for ( auto& fam : vFam ) {
        total_read_n += fam.read_n;

        for ( size_t t(0); t != fam.read_n; ++t ) {

            SeqLib::BamRecord rd;
            vecR[ fam.fh ].GetNextRecord(rd);

            if ( ok )
                w.WriteRecord(rd);
        }
    }

    if ( ok )
        outRdf << vFam[0].chr_start_end << "\t" << total_read_n << endl;
}

void add2mapRegionFamV( ifstream *vecF, ulong &i,  mRegFamV &rfmap, SeqLib::BamHeader &head )
{
    ifstream &fh = vecF[i];

    if ( !fh.good() ) {
        return;
    }

    string line, chr, start, end;
    Fam fam;
    fam.fh = i;

    if ( getline(fh,line) ) {

        istringstream stm(line);
        stm >> chr >> start >> end >> fam.read_n;
        fam.chr_start_end = chr + "\t" + start + "\t" + end;

        GenomicRegion rg(chr, start, end, head);
        rfmap[rg].push_back(fam);
    }
}

void mergeBlockFiles(SeqLib::BamWriter &w, ofstream &outRdf, string &outBam, ulong &part, BamHeader &head, ulong &rdfSizeCut )
{
    mRegFamV rfmap;
    vector<SeqLib::BamReader> vecR;
    ifstream vecF[part];
    string rmBlockFiles = "rm -f ";

    for ( size_t i(0); i != part; i++ ) {
        string outBam_p = outBam + "." + itoa(i);
        string outRdf_p = outBam_p + ".rdf";
        rmBlockFiles += outBam_p + " " + outRdf_p + " ";

        SeqLib::BamReader bam;
        bam.Open(outBam_p);
        vecR.push_back(bam);

        vecF[i].open( outRdf_p.c_str() );
        if ( !vecF[i].is_open() )  cerr << "cannot open file: " << outRdf_p << endl;

        add2mapRegionFamV(vecF, i, rfmap, head);
    }

    for (;;) {
        if ( rfmap.empty() ) break;

        mRegFamV::iterator itr = rfmap.begin();
        vector<Fam> &vFam = itr->second;
        ulong reads_num = sum_reads_n(vFam);

        if ( reads_num >= rdfSizeCut )
            show(w, outRdf, vecR, vecF, vFam, true);
        else
            show(w, outRdf, vecR, vecF, vFam, false);

        vector<ulong> id;
        for ( auto& fam : vFam )  id.push_back(fam.fh);

        rfmap.erase( itr );

        for ( auto& i : id ) add2mapRegionFamV(vecF, i, rfmap, head);
    }

    for ( size_t i(0); i != part; ++i ) {
        vecF[i].close();
        vecR[i].Close();
    }

    w.Close();
    outRdf.close();

    cerr << rmBlockFiles << endl;
    system( rmBlockFiles.c_str() );
}

void printBlockFile(string &outBam, ulong &part, mRegBrV &h, const BamHeader &head )
{
    string outBam_p = outBam + "." + itoa(part);
    string outRdf_p = outBam_p + ".rdf";

    cerr << "output part " << part << endl;

    SeqLib::BamWriter bam(SeqLib::BAM);
    if ( !bam.Open(outBam_p) )
        cerr << "cannot open file: " << outBam_p << endl;

    bam.SetHeader(head);
    bam.WriteHeader();

    ofstream rdf( outRdf_p.c_str() );
    if ( !rdf.is_open() )
        cerr << "cannot open file: " << outRdf_p << endl;

    for ( mRegBrV::iterator it = h.begin(); it != h.end(); it++ ) {
        uint32_t vecSize = it->second.size();

        const GenomicRegion &gr = it->first;

        rdf << gr.ChrName(head)  << "\t"
            << gr.pos1 + 1       << "\t"
            << gr.pos2 + 1       << '\t'
            << vecSize           << endl;

        for ( auto& rd : it->second )
            bam.WriteRecord(rd);
    }

    bam.Close();
    rdf.close();
}

