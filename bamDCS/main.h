/*
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#define PROGRAM "bamDCS"
#define VERSION "v2.4"
#define AUTHORS "yerui"
#define CONTACT "yerui@connect.hku.hk"
#define REMARKS "(generate double strand consensus reads for low-frequency mutation)"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "gzstream.h"
#include "boost/math/distributions/chi_squared.hpp"

using namespace __gnu_cxx;
using namespace std;
using namespace SeqLib;

typedef unsigned long ulong;

class Option {
    public:
        Option():
            baseQuaCutoff(20),
            mapQuaCutoff(30),
            minSupOnEachStrand(3),
            maxSupOnEachStrand(3000),
            Ncutoff(0.1),
            minFractionInFam(0.002),
            freqPrecision(0.00001),
            lhrGapCutoff(2.0),
            phredOffset(33),
            minSupOnHaplotype(3),
            filtSoftClip(false),
            outBamFile(""),
            sscsOut(false),
            singleOut(false),
            debug(false),
            pvalue(0.001),
            pcrError(1.0e-5),
            softEndTrim(5),
            randNread(0) {}

        ~Option(){}

        int baseQuaCutoff;
        int mapQuaCutoff;
        ulong minSupOnEachStrand;
        ulong maxSupOnEachStrand;

        double Ncutoff;
        double minFractionInFam;
        double freqPrecision;
        double lhrGapCutoff;

        int phredOffset;
        ulong minSupOnHaplotype;
        bool filtSoftClip;
        string outBamFile;
        bool sscsOut;
        bool singleOut;
        bool debug;
        double pvalue;
        double pcrError;
        ulong softEndTrim;
        size_t randNread;
};

typedef map<string, vector<SeqLib::BamRecord> > mStrBrV;
typedef map<char,double> mCharDouble;
typedef map<char, ulong> mCharUlong;
typedef map<string,ulong> mStrUlong;
typedef map<char, vector<double> > mCvD;

typedef pair<char,ulong> pCharUlong;
typedef pair<double, set<char> > pDoubleCharSet;

typedef vector< set<char> > vCharSet;
typedef vector<string> vString;

// descending sort pair
inline bool _cmpByFirst(const pDoubleCharSet &a, const pDoubleCharSet &b) { return a.first > b.first; }
inline bool _cmpBySecond(const pCharUlong &a, const pCharUlong &b) { return a.second > b.second; }

inline double errorRate(const char &q, const Option &opt);
inline vector<double> quaToErrorRate(const string &qs, const Option &opt);
inline string _itoa( size_t &i );
inline ulong countN(const string &s);

inline bool lowQuality(const char &, const Option &);
inline string reverseComplement(const string &);

void usage();
SeqLib::BamHeader removeReadGroup( const SeqLib::BamHeader & );
void calibrateFam(mStrBrV &);
bool testCigarFam(mStrBrV &, mStrBrV &, const string &, const Option &);
void printConsensusRead(ogzstream &,ogzstream &,mStrBrV &,mStrBrV &,const string &,const Option &,const string &, SeqLib::BamWriter &);
void sscs(mStrBrV &,mStrBrV &,const string &,const Option &,const string &, SeqLib::BamWriter &);
void addPcrError( vector< mCvD > &wcHetPos, const Option &opt);

string getQuaFromPvalue( const vector< mCvD > &quaV, const string &s, const Option &opt );
string adjust_p(const string &qs, const Option &opt);

vector< mCvD > hetPoint(const BamRecordVector &, const Option &);

vector< mCvD > zipHetPoint(const vector< mCvD > &, const vector< mCvD > &, const Option &opt);
map<string, long> consensusSeq(const BamRecordVector &, const BamRecordVector &, const vector< mCvD > &, const Option &);
map<string, long> consensusSeq(const BamRecordVector &wcBrV, const vector< mCvD > &wcHetPos, const Option &opt);

mCvD  llh_genotype(const string &, const string &, const Option &);

string ignoreError(const string &, const vector< mCvD > &);
string ignoreError2(const string &s, const vector< mCvD > &v);

void trimEnd(SeqLib::BamRecord &br, const Option &opt);
size_t getGenomePosition(size_t pos, ulong i, const Cigar &cg);

// sub functions

void usage() {
    cout << "Program: " << PROGRAM << "  " << REMARKS << "\n"
        "Version: " << VERSION << "\n"
        "Authors: " << AUTHORS << "\n"
        "Contact: " << CONTACT << "\n\n"

        "Options: " << PROGRAM << " in.bam out_prefix\n\n"

        "    -q [i]     base quality cutoff [20]\n"
        "    -Q [i]     map quality cutoff [30]\n"

        "    -s [i]     min support num on each strand [3]\n"
        "    -S [i]     max support num on each strand [3000]\n"

        "    -N [f]     max fraction of N on consensus read [0.1]\n"
        "    -f [f]     min fraction of alterative allele in a read family [0.002]\n"

        "    -e [f]     precision of allele frequency while calculate likelihood ratio [0.00001]\n"
        "    -g [f]     gap between likelihood ratios of major and rest genotypes [2.0]\n"
        "    -x [i]     Encoding offset for phred quality scores [33]\n"
        "    -t [i]     min support num to construct a haplotype seq [3]\n"
        "    -c         discard soft-clipping reads [false]\n"
        "    -C [i]     soft trim N base from both ends of read [5]\n"
        "    -p [f]     expected PCR error rate during library construction [1e-5]\n"
        "    -n [i]     N read pairs to be randomly selected and used [0, total]\n"

        "    -o [s]     output bam File directly []\n"
        "    -a         output single pair-end reads [false]\n"
        "    -b         output SSCS [false]\n"
        "    -d         debug mode [false]\n"
        "    -h         help\n"
        "    -v         version\n"
        "\n";
}

SeqLib::BamHeader removeReadGroup( const SeqLib::BamHeader &oldHead )
{
    istringstream itm( oldHead.AsString() );
    string line;
    string s("");

    while (getline(itm, line) ) {
        size_t i = line.find("@RG");

        if ( i == 0 ) {
            s += "@RG\tID:foo\tSM:bar\n";
        }
        if ( i != 0 ) {
            s += line + "\n";
        }
    }

    // add rg if the old head does not have rg
    if ( s.find("@RG") == string::npos ) {
        s += "@RG\tID:foo\tSM:bar\n";
    }

    SeqLib::BamHeader newHead(s);
    return newHead;
}

inline double errorRate(const char &q, const Option &opt)
{
    return pow( 10, (double)(opt.phredOffset-(int)q)/10 );
}

inline vector<double> quaToErrorRate(const string &qs, const Option &opt)
{
    vector<double> eV;
    for ( size_t i(0); i != qs.size(); i++ )  eV.push_back( errorRate(qs[i], opt) );
    return eV;
}

inline string _itoa( size_t &i )
{
    ostringstream osm;
    osm << i;
    return osm.str();
}

inline ulong countN(const string &s)
{
    ulong n(0);
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) n++;
    }
    return n;
}

inline bool lowQuality(const char &q, const Option &opt)
{
    return ( int(q) - opt.phredOffset < opt.baseQuaCutoff );
}

inline string reverseComplement(const string &seq)
{
    string str("");
    for ( string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); it++) {
        switch (*it) {
            case 'A': str += "T";    break;
            case 'C': str += "G";    break;
            case 'G': str += "C";    break;
            case 'T': str += "A";    break;
            case 'N':
            default:  str += "N";    break;
        }
    }

    return str;
}

inline char errorRateToChar(const double &q, const Option &opt)
{
    return char( opt.phredOffset - (int)(-10 * log(q) / log(10.0) ) );
}

void calibrateFam(mStrBrV &)
{
    return;
}

// check family size on cigar
bool testCigarFam(mStrBrV &watsonFam, mStrBrV &crickFam, const string &cg, const Option &opt)
{
    mStrBrV::const_iterator w = watsonFam.find( cg );
    mStrBrV::const_iterator c =  crickFam.find( cg );

    if (    w != watsonFam.end()
            && c !=  crickFam.end()
            && w->second.size() >= opt.minSupOnEachStrand * 2
            && c->second.size() >= opt.minSupOnEachStrand * 2
            && w->second.size() <= opt.maxSupOnEachStrand * 2
            && c->second.size() <= opt.maxSupOnEachStrand * 2
       )
        return true;
    else
        return false;
}

// output consensus Read pair into fq files
void printConsensusRead(
        ogzstream & fq1,
        ogzstream & fq2,
        mStrBrV & watsonFam,
        mStrBrV & crickFam,
        const string & cg,
        const Option & opt,
        const string & chrBegEnd,
        SeqLib::BamWriter & writer
        ) {
    // length of read1 & read2 are same, So connect them to simplify workflow

    // 0-based index of heterozygous point --> vector of allele set
    vector< mCvD > wHetPos, cHetPos, sameHetPos;

    // find heterzygous point on each fam
    wHetPos = hetPoint(watsonFam[cg], opt);
    cHetPos = hetPoint( crickFam[cg], opt);

    // find consistent het point by comparing watson & crick family
    sameHetPos = zipHetPoint(wHetPos, cHetPos, opt);

    // watson & crick should be concordant on hom point. if not, set N
    map<string, long> mSeqN = consensusSeq(watsonFam[cg], crickFam[cg], sameHetPos, opt);

    string Qname("@");
    Qname += chrBegEnd + ":" + cg + ":";

    if ( mSeqN.empty() ) return;

    double dcsN(0.0);
    for ( auto p : mSeqN ) dcsN += (double)p.second;

    int i(0);
    for ( auto p : mSeqN ) {
        i++;
        string seq = p.first;
        double frac( p.second/dcsN );

        size_t length = seq.size()/2;
        string quaStr = getQuaFromPvalue( sameHetPos, seq, opt );

        string rd1 = seq.substr( 0, length );
        string rd2 = reverseComplement( seq.substr(length) );

        string quaStr1 = quaStr.substr( 0, length );
        string quaStr2 = quaStr.substr( length    );

        reverse( quaStr2.begin(), quaStr2.end() );

        if ( opt.outBamFile.size() > 0 ) {
            SeqLib::BamRecord br1 = watsonFam[cg].at(0);
            SeqLib::BamRecord br2 = watsonFam[cg].at(1);

            if ( opt.softEndTrim > 0 ) {
                trimEnd(br1, opt);
                trimEnd(br2, opt);
            }

            br1.SetSequence( rd1 );
            br2.SetSequence( seq.substr(length) );

            if ( br1.Position() < 1 || br2.Position() < 1 )
                continue;

            reverse( quaStr2.begin(), quaStr2.end() );
            br1.SetQualities( quaStr1, 33 );
            br2.SetQualities( quaStr2, 33 );

            br1.RemoveAllTags();
            br2.RemoveAllTags();

            br1.AddZTag("RG", "foo");
            br2.AddZTag("RG", "foo");

            if ( frac > 0 ) {
                br1.AddZTag("fr", to_string(frac));
                br2.AddZTag("fr", to_string(frac));

                ostringstream famSize;
                famSize << watsonFam[cg].size() << "," << crickFam[cg].size();
                br1.AddZTag("fs", famSize.str() );
                br2.AddZTag("fs", famSize.str() );
            }

            if ( mSeqN.size() > 1 ) {
                br1.AddZTag("sp", to_string(i));
                br2.AddZTag("sp", to_string(i));
            }

			if ( br1.GetCigar().NumQueryConsumed() != br1.Length() ) {
				cerr << "WARN: cigar and seq length differ for " << br1 << endl;
				continue;
			}

			if ( br2.GetCigar().NumQueryConsumed() != br2.Length() ) {
				cerr << "WARN: cigar and seq length differ for " << br2 << endl;
				continue;
			}

			try {
				writer.WriteRecord( br1 );
				writer.WriteRecord( br2 );
			}
			catch ( std::runtime_error &e ) {
				cerr << "WARN: WriteRecord error, skip\n" << br1 << "\n" << br2 << endl;
			}
        }
    }
}

void sscs(mStrBrV &watsonFam,
        mStrBrV &crickFam,
        const string &cg,
        const Option &opt,
        const string &chrBegEnd,
        SeqLib::BamWriter &writer )
{
    vector< mCvD > wcHetPos;
    BamRecordVector wcBrV;
    size_t w_size(0), c_size(0);

    // wcBrV contain br from both Fam
    for ( auto &br : watsonFam[cg] ) {
        wcBrV.push_back(br);
        w_size++;
    }

    for ( auto &br : crickFam[cg] ) {
        wcBrV.push_back(br);
        c_size++;
    }

    // only one pair of reads, output directly
    if ( opt.singleOut && wcBrV.size() == 2 && opt.outBamFile.size() > 0 ) {
        if ( opt.softEndTrim > 0 ) {
            trimEnd( wcBrV[0], opt );
            trimEnd( wcBrV[1], opt );
        }

		ostringstream famSize;
		famSize << w_size << "," << c_size;
		wcBrV[0].AddZTag("fs", famSize.str() );
		wcBrV[1].AddZTag("fs", famSize.str() );

		writer.WriteRecord( wcBrV[0] );
		writer.WriteRecord( wcBrV[1] );
		return;
	}
/*
	// skip sscs if there is indel and not supported by both strands
	size_t insIt = cg.find('I');
    size_t delIt = cg.find('D');
    if (   ( insIt != string::npos or delIt != string::npos )
        && ( w_size == 0 or c_size == 0 )
       ) {
        return;
    }
*/
    // more than one pairs of reads, calculate pvalue + PCR error
    wcHetPos = hetPoint(wcBrV, opt);
    addPcrError(wcHetPos, opt);

    map<string, long> mSeqN = consensusSeq(wcBrV, wcHetPos, opt);

    string Qname("@");
    Qname += chrBegEnd + ":" + cg + ":";

    if ( mSeqN.empty() ) return;

    double dcsN(0);
    for ( auto p : mSeqN ) dcsN += (double)p.second;

	int i(0);
    for ( auto p : mSeqN ) {
		i++;
        string seq = p.first;
        double frac( p.second/dcsN );

        size_t length = seq.size()/2;
        string quaStr = getQuaFromPvalue( wcHetPos, seq, opt );

        string rd1 = seq.substr( 0, length );
        string rd2 = reverseComplement( seq.substr(length) );

        string quaStr1 = quaStr.substr( 0, length );
        string quaStr2 = quaStr.substr( length    );

        reverse( quaStr2.begin(), quaStr2.end() );
        if ( opt.outBamFile.size() > 0 ) {
            SeqLib::BamRecord br1 = wcBrV.at(0);
            SeqLib::BamRecord br2 = wcBrV.at(1);

            if ( opt.softEndTrim > 0 ) {
                trimEnd(br1, opt);
                trimEnd(br2, opt);
            }

            br1.SetSequence( rd1 );
            br2.SetSequence( seq.substr(length) );

            reverse( quaStr2.begin(), quaStr2.end() );
            br1.SetQualities( quaStr1, 33 );
            br2.SetQualities( quaStr2, 33 );

            br1.RemoveAllTags();
            br2.RemoveAllTags();

            br1.AddZTag("RG", "foo");
            br2.AddZTag("RG", "foo");

            if (frac > 0) {
                br1.AddZTag("fr", to_string(frac));
                br2.AddZTag("fr", to_string(frac));

                ostringstream famSize;
                famSize << w_size << "," << c_size;
                br1.AddZTag("fs", famSize.str() );
                br2.AddZTag("fs", famSize.str() );

            }

            if ( mSeqN.size() > 1 ) {
                br1.AddZTag("sp", to_string(i));
                br2.AddZTag("sp", to_string(i));
            }

			if ( br1.GetCigar().NumQueryConsumed() != br1.Length() ) {
				cerr << "WARN: cigar and seq length differ for " << br1 << endl;
				continue;
			}

			if ( br2.GetCigar().NumQueryConsumed() != br2.Length() ) {
				cerr << "WARN: cigar and seq length differ for " << br2 << endl;
				continue;
			}


            writer.WriteRecord( br1 );
            writer.WriteRecord( br2 );
        }
    }

    return;
}

// get heterozygous points based on watson or crick family only
vector< mCvD > hetPoint(const BamRecordVector &brV, const Option &opt)
{
    vector< mCvD > pt;
    vString seqV, quaV;

    // too few reads to support heterozygous point. two reads form a pair.
    if ( brV.size() < opt.minSupOnEachStrand * 2 )
        return pt;

    string seq(""), qua("");
    for ( auto & br : brV ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
            qua = br.Qualities();
        }
        else { // - strand
            seq += br.Sequence();
            qua += br.Qualities();

            seqV.push_back(seq);
            quaV.push_back(qua);
        }
    }

    size_t length = seq.size();

    // fetch each column of alleles and quality scores
    for ( size_t j(0); j != length; j++ ) {
        string base(""), qual("");

        for ( size_t i(0); i != seqV.size(); i++ ) {
            base += seqV[i][j];
            qual += quaV[i][j];
        }

        pt.push_back( llh_genotype(base, qual, opt) );
    }

    return pt;
}

void addPcrError( vector< mCvD > &wcHetPos, const Option &opt)
{
    for ( auto &ntPF : wcHetPos ) {
        for ( mCvD::iterator it = ntPF.begin(); it != ntPF.end(); it++ ) {
            it->second.at(0) += 10 * opt.pcrError;
            if ( it->second.at(0) > 1 ) it->second.at(0) = 1.0;
        }
    }

    return;
}

string getQuaFromPvalue( const vector< mCvD > &quaV, const string &s, const Option &opt )
{
    ostringstream q("");

    for ( size_t i(0); i != s.size(); i++ ) {

        if ( s[i] == 'N' ) {
            q << '$';
        }
        else {
            mCvD::const_iterator it = quaV[i].find( s[i] );

            if ( it != quaV[i].end() ) {
                double f = -10.0 * log( it->second.at(0) )/log(10.0);
                q << (char)( opt.phredOffset + int(f) );
            }
            else {
                q << '$';
            }

        }
    }

    return q.str();
}

vector< mCvD > zipHetPoint(const vector< mCvD > &w, const vector< mCvD > &c, const Option &opt)
{
    vector< mCvD > samePt;

    if ( w.empty() || c.empty() )
        return samePt;

    for ( size_t j(0); j != w.size(); j++ ) {
        mCvD nt;

        for ( mCvD::const_iterator wi = w[j].begin(); wi != w[j].end(); wi++ ) {
            if ( wi->second.at(0) > opt.pvalue ) continue;

            mCvD::const_iterator ci = c[j].find( wi->first );
            if ( ci != c[j].end() ) {
                if ( ci->second.at(0) > opt.pvalue ) continue;

                nt[ wi->first ].push_back( wi->second.at(0) + ci->second.at(0) - wi->second.at(0) * ci->second.at(0) + 10 * pow(opt.pcrError,2) );

                nt[ wi->first ].push_back( (wi->second.at(1) + ci->second.at(1)) / 2);
            }
        }

        samePt.push_back(nt);
    }

    return samePt;
}

map<string, long> consensusSeq(const BamRecordVector &w, const BamRecordVector &c, const vector< mCvD > &sameHetPos, const Option &opt)
{
    map<string, long> mSeqN;
    if ( sameHetPos.empty() ) return mSeqN;

    mStrUlong mSeqN_w, mSeqN_c;
    string seq("");

    // count N
    size_t Ncnt(0);
    for ( auto & p : sameHetPos ) {
        if ( p.empty() )
            Ncnt++;
    }

    if ( Ncnt > sameHetPos.size() * opt.Ncutoff )
        return mSeqN;

    // get potential seq based on original read and heterozygous site
    for ( auto &br : w ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
        }
        else { // - strand
            seq += br.Sequence();
            mSeqN_w[ ignoreError(seq, sameHetPos) ]++;
        }
    }

    for ( auto &br : c ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
        }
        else { // - strand
            seq += br.Sequence();
            mSeqN_c[ ignoreError(seq, sameHetPos) ]++;
        }
    }

    for ( mStrUlong::iterator it = mSeqN_w.begin(); it != mSeqN_w.end(); it++ ) {
        if ( it->second >= opt.minSupOnHaplotype ) {
            size_t len = it->first.size() / 2;
            if (   countN( it->first.substr(  0,len) ) > opt.Ncutoff * len
                    && countN( it->first.substr(len,len) ) > opt.Ncutoff * len
               )
                continue;

            mStrUlong::iterator ct = mSeqN_c.find( it->first );

            if (   ct != mSeqN_c.end() && ct->second >= opt.minSupOnHaplotype )
                mSeqN[ it->first ] = it->second + ct->second;
        }
    }

    return mSeqN;
}

// for sscs
map<string, long> consensusSeq(const BamRecordVector &wcBrV,
        const vector< mCvD > &wcHetPos,
        const Option &opt )
{
    map<string, long> mSeqN;
    if ( wcHetPos.empty() ) return mSeqN;

    mStrUlong mSeqN_wc;
    string seq("");

    // count N
    size_t Ncnt(0);
    for ( auto & p : wcHetPos ) {
        if ( p.empty() )
            Ncnt++;
    }

    if ( Ncnt > wcHetPos.size() * opt.Ncutoff )
        return mSeqN;

    // get potential seq based on original read and heterozygous site
    for ( auto &br : wcBrV ) {
        if ( !br.ReverseFlag() ) { // + strand
            seq = br.Sequence();
        }
        else { // - strand
            seq += br.Sequence();
            mSeqN_wc[ ignoreError(seq, wcHetPos) ]++;
        }
    }

    for ( mStrUlong::iterator it = mSeqN_wc.begin(); it != mSeqN_wc.end(); it++ ) {
        if ( it->second >= opt.minSupOnHaplotype ) {
            size_t len = it->first.size() / 2;
            if (   countN( it->first.substr(0,len) ) > opt.Ncutoff * len
                    && countN( it->first.substr(len,len) ) > opt.Ncutoff * len
               )  continue;

            mSeqN[ it->first ] = it->second;
        }
    }

    return mSeqN;
}

string ignoreError(const string &s, const vector< mCvD > &v)
{
    string seq("");
    for ( size_t i(0); i != s.size(); i++ ) {

        if ( v[i].empty() ) {
            seq += "N";
        }
        else if ( v[i].size() == 1 ) {
            seq += v[i].begin()->first;
        }
        else {

            mCvD::const_iterator it = v[i].find( s[i] );
            if ( it != v[i].end() ) {
                seq += s[i];
            }
            else {
                seq += "N";
            }
        }
    }
    return seq;
}

string ignoreError2(const string &s, const vector< mCvD > &v)
{
    string seq("");
    for ( size_t i(0); i != s.size(); i++ ) {

        if ( v[i].empty() ) {
            seq += "N";
        }
        else if ( v[i].size() == 1 ) {
            seq += v[i].begin()->first;
        }
        else {
            seq += "N";
        }
    }
    return seq;
}

string adjust_p(const string &qs, const Option &opt)
{
    ostringstream o;

    map<double, vector<char> > m;
    vector<double> v = quaToErrorRate(qs, opt);

    if ( v.size() > 1 ) {
        sort( v.begin(), v.end() );

        ulong total( v.size() );
        for ( size_t i(0); i != total; i++ ) {
            m[ v[i] ].push_back( errorRateToChar(v[i], opt) );
        }

        for ( auto & q : qs) {
            double e = errorRate(q, opt);
            o << m[ e ].back();
            m[ e ].pop_back();
        }

        return o.str();
    }
    else {
        return qs;
    }
}

void trimEnd(SeqLib::BamRecord &br, const Option &opt)
{
    CigarField sc('S', opt.softEndTrim );
    Cigar cg = br.GetCigar();

    if ( cg.size() == 0 ) {
        return;
    }

    size_t head(0), tail(0);
    string ty("=XMIS");

    //trim head
    Cigar nc;

    //cerr << "old: " << cg << ' ' << br.Position() << ' ';

    for (size_t i(0); i != cg.size(); i++ ) {
        if ( ty.find( cg[i].Type() ) != string::npos ) {
            head += cg[i].Length();
        }

        if ( head >= opt.softEndTrim ) {
            nc.add(sc);
            size_t remain = head - opt.softEndTrim;
            if ( remain > 0 ) {
                CigarField tmp(cg[i].Type(), remain);
                nc.add(tmp);

                for (size_t j(i+1); j != cg.size(); j++) nc.add(cg[j]);
            }
            else {
                for (size_t j(i+1); j != cg.size(); j++) nc.add(cg[j]);
            }

            break;
        }
    }

    //trim tail
    Cigar rc;

    for ( int i( nc.size()-1 ); i >= 0; i-- ) {
        if ( ty.find( nc[i].Type() ) != string::npos ) {
            tail += nc[i].Length();
        }

        if ( tail >= opt.softEndTrim ) {
            rc.add(sc);

            int remain = tail - opt.softEndTrim;
            if ( remain > 0 ) {
                CigarField tmp(nc[i].Type(), remain);
                rc.add(tmp);
                for (int j(i-1); j >= 0; j--) rc.add(nc[j]);
            }
            else {
                for (int j(i-1); j >= 0; j--) rc.add(nc[j]);
            }

            break;
        }
    }

    //cerr << " new: " << rc << ' ' << br.Position() << ' ';

    // deal with I D M near head & tail

    vector<CigarField> vrc;
    if ( rc[1].Type() == 'I' || rc[1].Type() == 'S' ) { // combine I & S
        CigarField tf('S', rc[0].Length() + rc[1].Length() );
        vrc.push_back(tf);
        for ( size_t i(2); i < rc.size(); i++ ) vrc.push_back( rc[i] );
    }
    else if ( rc[1].Type() == 'D' || rc[2].Type() == 'N' ) { // skip D
        vrc.push_back( rc[0] );
        for ( size_t i(2); i < rc.size(); i++ ) vrc.push_back( rc[i] );
    }
    else {
        for ( auto &p : rc ) vrc.push_back(p);
    }

    reverse( vrc.begin(), vrc.end() );

    // reverse vrc
    vector<CigarField> rv;
    if ( vrc[1].Type() == 'I' || rc[1].Type() == 'S' ) { // combine I & S
        CigarField tf('S', vrc[0].Length() + vrc[1].Length() );
        rv.push_back(tf);
        for ( size_t i(2); i < vrc.size(); i++ ) rv.push_back( vrc[i] );
    }
    else if ( vrc[1].Type() == 'D' || vrc[2].Type() == 'N' ) { // skip D

        rv.push_back( vrc[0] );
        for ( size_t i(2); i < vrc.size(); i++ ) rv.push_back( vrc[i] );
    }
    else {
        for ( auto &p : vrc ) rv.push_back(p);
    }

    // reverse
    Cigar rev;

    for ( auto & p : rv ) {
        rev.add( p );
    }

    // set new start position of reads
    size_t rev_S_length = ( rev[0].Type() == 'S' ? rev[0].Length() : 0 );

    br.SetPosition( getGenomePosition(br.Position(), rev_S_length, cg) );

    //cerr << " final: " << rev << ' ' << br.Position() << endl;

    br.SetCigar( rev );
}

size_t getGenomePosition(size_t pos, ulong i, const Cigar &cg)
{
    for ( vector<CigarField>::const_iterator it = cg.begin(); it != cg.end(); it++ ) { 
        char t = it->Type();
        size_t n = it->Length();

        switch (t) {
            case '=':
            case 'X':
            case 'M':
                if (i < n)  return (pos+i);
                else        (pos+=n, i-=n);
                break;
            case 'I':
                if (i < n)  return pos;
                else        (i -= n); 
                break;
            case 'S':
                if (i < n)  return -1; 
                else        i -= n;
                break;

            case 'N':
            case 'D':       pos += n;  break;
            case 'H':
            case 'P':                  break;

            default:  cerr << "unknown cigar field: " << cg << endl, exit(1);
        }
    }

    return pos + i;
}
#include <boost/math/distributions/chi_squared.hpp>

// needed by alglib
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "ap.h"

struct fn_data {
    std::string         base;
    std::vector<double> errRateV;
};

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
    double l_c(0.0);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i]/3;
        switch( base[i] ) {
            case 'A': l_c += log( (1-4*e) * theta[0] + e ); break;
            case 'C': l_c += log( (1-4*e) * theta[1] + e ); break;
            case 'G': l_c += log( (1-4*e) * theta[2] + e ); break;
            case 'T': l_c += log( (1-4*e) * theta[3] + e ); break;
            default: cerr << "unknown base in " << base << endl,exit(1);
        }
    }

    return l_c;
}

// -composite score function: -U_c(theta)
// return a column vector
//
alglib::real_1d_array composite_score (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
    alglib::real_1d_array U_c = "[0,0,0,0]";

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i]/3;

        switch( base[i] ) {
            case 'A': U_c[0] -= (1-4*e) / ( (1-4*e)*theta[0] + e );  break;
            case 'C': U_c[1] -= (1-4*e) / ( (1-4*e)*theta[1] + e );  break;
            case 'G': U_c[2] -= (1-4*e) / ( (1-4*e)*theta[2] + e );  break;
            case 'T': U_c[3] -= (1-4*e) / ( (1-4*e)*theta[3] + e );  break;
            default: cerr << "unknown base in " << base << endl, exit(1);
        }
    }

    return U_c;
}

// gradient optimization
void function1_grad (
        const alglib::real_1d_array  &x,
        double                       &func,
        alglib::real_1d_array        &grad,
        void                         *opt_data )
{
    fn_data* objfn_data = reinterpret_cast<fn_data*>(opt_data);
    const std::string &base = objfn_data->base;
    const std::vector<double> &errRateV = objfn_data->errRateV;

    func = -composite_LogLikelihood( base, errRateV, x);
    grad = composite_score( base, errRateV, x );
}

string initAlleleFreq (
        mCharUlong       &fr,
        const double     &depth,
        const char       &except_b )
{
    ostringstream s;
    s << '[';

    if ( depth == fr[except_b] ) {
        for ( auto b : "ACGT" ) {
            b == except_b ? s << 0.0 : s << 0.333333333;
            b == 'T' ? s << ']' : s << ',';
            if ( b == 'T' ) break;
        }
    }
    else {
        for ( auto b : "ACGT" ) {
            b == except_b ? s << 0.0 : s << fr[b] / (depth - fr[except_b]);
            b == 'T' ? s << ']' : s << ',';
            if ( b== 'T' ) break;
        }
    }
    return s.str();
}

string initAlleleFreq (
        mCharUlong       &fr,
        double           depth,
        const set<char>  &except_bs )
{
    ostringstream s;

    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        if ( it != except_bs.end() ) {
            fr[b] = 0;
            depth -= fr[b];
        }

        if ( b== 'T' ) break;
    }

    if ( depth == 0 ) {
        s  << "[0.25,0.25,0.25,0.25]";
        return s.str();
    }

    s << '[';
    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        it != except_bs.end() ? s << 0.0 : s << fr[b] / depth;  // depth here don't contain num of except_bs
        b == 'T' ? s << ']' : s << ',';
        if ( b== 'T' ) break;
    }

    return s.str();
}

string _upBoundary(const char except_b)
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        b == except_b ? s << 0 : s << 1;
        b == 'T' ? s << ']' : s << ',';
        if ( b == 'T' ) break;
    }
    return s.str();
}

string _upBoundary(const set<char> except_bs)
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        set<char>::const_iterator it = except_bs.find(b);
        it != except_bs.end() ? s << 0 : s << 1;
        b == 'T' ? s << ']' : s << ',';
        if ( b == 'T' ) break;
    }
    return s.str();
}

map<char, vector<double> > llh_genotype(const string &s, const string &q, const Option &opt)
{
    map<char, vector<double> > ntPF; // nt => pvalue, fraction
    boost::math::chi_squared X2_dist(1);

    mCharUlong fr;
    for ( auto b : "ACGTN" ) {
        fr[b] = 0;
        if ( b == 'N' ) break;
    }

    string new_s(""), new_q("");
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' || s[i] == '*' )  continue; // fr['N'] == 0
        fr[ s[i] ]++;
        new_s += s[i];
        new_q += q[i];
    }

    double depth(new_s.size());
    vector<double> errV = quaToErrorRate(new_q, opt);

    if ( depth == 0 ) return ntPF;

    // quick return if equal or less then one allele
    if ( fr.empty() ) {
        return ntPF;
    }

    // sort by frequent
    vector<pair<char, ulong> > ntV( fr.begin(), fr.end() );

    // only one allele
    if ( ntV.size() == 1 ) {
        if ( ntV[0].second >= opt.minSupOnEachStrand ) {
            ntPF[ ntV[0].first ].push_back(0.0);
            ntPF[ ntV[0].first ].push_back(1.0);
        }

        return ntPF;
    }

    if ( ntV.size() > 1 )
        sort( ntV.begin(), ntV.end(), _cmpBySecond ); // descending sort

    // delete pair<allele, supportNum> which has too few support reads or too small fraction
    while ( ntV.size() ) {
        if ( ntV.back().second < opt.minSupOnEachStrand || ntV.back().second / depth < opt.minFractionInFam )
            ntV.pop_back();
        else
            break;
    }

    // none allele remain
    if ( ntV.empty() ) {
        return ntPF;
    }

    // only one allele
    if ( ntV.size() == 1 ) {
        if ( ntV[0].second >= opt.minSupOnEachStrand ) {
            ntPF[ ntV[0].first ].push_back(0.0);
            ntPF[ ntV[0].first ].push_back(1.0);
        }

        return ntPF;
    }



    // more than one allele
    fn_data data;
    data.base = new_s;
    data.errRateV = errV;

    // var for alglib
    alglib::minbleicstate state;
    alglib::minbleicreport rep;
    double epsg(0.0001);
    double epsf(0.0);
    double epsx(0.0);
    alglib::ae_int_t maxits(0);

    // constraint: sum of frequency of 4 alleles == 1
    alglib::real_2d_array c = "[[1,1,1,1,1]]";
    alglib::integer_1d_array ct = "[0]";
    alglib::real_1d_array bndl = "[0,0,0,0]";

    // four allele maximize
    double cl_4(0.0);
    try {
        string AFstr = initAlleleFreq(fr, depth, 'N');
        alglib::real_1d_array alg_x = AFstr.c_str();
        alglib::real_1d_array bndu = "[1,1,1,1]";

        alglib::minbleiccreate(alg_x, state);
        alglib::minbleicsetlc(state, c, ct);
        alglib::minbleicsetbc(state, bndl, bndu);
        alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
        alglib::minbleicoptimize(state, function1_grad, NULL, &data );
        alglib::minbleicresults(state, alg_x, rep);

        if ( opt.debug ) {
            printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
            printf("%s\n", alg_x.tostring(20).c_str());
        }

        cl_4 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
        if ( opt.debug ) cout << "cl_4: " << setprecision(20) << cl_4 << endl;
    }
    catch ( alglib::ap_error &e ) {
        cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q << "]" << endl;
    }

    map<char, string> init_V;
    map<char, string> bndu_V;

    for ( auto b : "ACGT" ) {
        init_V[b] = initAlleleFreq(fr, depth, b);
        bndu_V[b] = _upBoundary(b);
        if ( b == 'T' ) break;
    }

    for ( mCharUlong::const_iterator it = fr.begin(); it != fr.end(); it++ )
    {
        if ( it->second < opt.minSupOnEachStrand || it->second/depth < opt.minFractionInFam ) continue;

        double cl_3(0.0);
        try {
            alglib::real_1d_array alg_x = init_V[ it->first ].c_str();
            alglib::real_1d_array bndu = bndu_V[ it->first ].c_str();

            alglib::minbleiccreate(alg_x, state);
            alglib::minbleicsetlc(state, c, ct);
            alglib::minbleicsetbc(state, bndl, bndu);
            alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
            alglib::minbleicoptimize(state, function1_grad, NULL, &data );
            alglib::minbleicresults(state, alg_x, rep);

            if ( opt.debug ) {
                printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
                printf("%s\n", alg_x.tostring(20).c_str());
            }

            cl_3 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
            if ( opt.debug ) cout << "cl_3: " << cl_3 << endl;
        }
        catch ( alglib::ap_error &e ) {
            cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q
                << "] for base[" << it->first << "]" << endl;
        }

        if ( cl_4 - cl_3 > opt.lhrGapCutoff )
            ntPF[ it->first ].push_back( 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3)) );
    }

    if ( ntPF.size() == 1 ) {
        ntPF[ ntPF.begin()->first ].push_back(1.0);
        return ntPF;
    }
    else if ( ntPF.size() > 1 ) {
        set<char> except_bs;
        for ( auto b : "ACGT" ) {
            map<char, vector<double> >::const_iterator it = ntPF.find(b);
            if ( it == ntPF.end() ) {  // not in ntPF
                except_bs.insert(b);
            }
            if ( b == 'T' ) break;
        }

        string AFstr = initAlleleFreq(fr, depth, except_bs);
        alglib::real_1d_array alg_x = AFstr.c_str();
/*
        try {
            string upBnd = _upBoundary(except_bs);
            alglib::real_1d_array bndu = upBnd.c_str();

            alglib::minbleiccreate(alg_x, state);
            alglib::minbleicsetlc(state, c, ct);
            alglib::minbleicsetbc(state, bndl, bndu);
            alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
            alglib::minbleicoptimize(state, function1_grad, NULL, &data );
            alglib::minbleicresults(state, alg_x, rep);

            if ( opt.debug ) {
                printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
                printf("%s\n", alg_x.tostring(20).c_str());
            }

            cl_4 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
            if ( opt.debug ) cout << "cl_4: " << setprecision(20) << cl_4 << endl;
        }
        catch ( alglib::ap_error &e ) {
            cerr << "catch error: " << e.msg << " at seq[" << new_s << "] qua[" << new_q << "]" << endl;
        }
*/
        string st = "ACGT";
        map<char, double> mBaseFrac;
        for ( int i(0); i != 4; i++ ) {
            mBaseFrac[ st[i] ] = alg_x[i];
        }

        for ( auto &p : ntPF ) {
            ntPF[ p.first ].push_back( mBaseFrac[p.first] );
        }

        return ntPF;
    }
    else if ( ntPF.size() > 4 ) {
        cerr << "ntPF contain unknown base" << endl, exit(1);
    }
    else {
        return ntPF;
    }
}

#endif // MAIN_H_
