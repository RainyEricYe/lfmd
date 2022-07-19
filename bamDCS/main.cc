/*
 */

#include "main.h"

int main( int argc, char **argv )
{

    if ( argc < 2 ) usage(), exit(1);

    srand (time(NULL));
    int c;
    Option opt;
    while ( (c=getopt(argc,argv,"q:Q:s:S:N:f:e:g:x:t:o:C:p:n:cabdvh")) != -1 ) {
        switch (c) {
            case 'q': opt.baseQuaCutoff = atoi(optarg);               break;
            case 'Q': opt.mapQuaCutoff  = atoi(optarg);               break;

            case 's': opt.minSupOnEachStrand = atoi(optarg);          break;
            case 'S': opt.maxSupOnEachStrand = atoi(optarg);          break;

            case 'N': opt.Ncutoff = atof(optarg);                     break;
            case 'f': opt.minFractionInFam      = atof(optarg);       break;

            case 'e': opt.freqPrecision = atof(optarg);               break;
            case 'g': opt.lhrGapCutoff  = atof(optarg);               break;
            case 'x': opt.phredOffset   = atoi(optarg);               break;

            case 't': opt.minSupOnHaplotype = atoi(optarg);           break;
            case 'o': opt.outBamFile = optarg;                        break;
            case 'c': opt.filtSoftClip = true;                        break;
            case 'C': opt.softEndTrim = atoi(optarg);                 break;
            case 'p': opt.pcrError = atof(optarg);                    break;
            case 'n': opt.randNread = atoi(optarg);                   break;

            case 'a': opt.singleOut = true;                           break;
            case 'b': opt.sscsOut = true;                             break;
            case 'd': opt.debug = true;                               break;
            case 'v': cerr << VERSION << endl;                        exit(1);
            case 'h':
            default:  usage();                                        exit(1);
        }
    }

    if ( opt.singleOut ) opt.sscsOut = true;

    if ( argc < optind + 2 ) usage(), exit(1);

    string  inBamFile( argv[optind] );
    string  inRdfFile = inBamFile + ".rdf";

    string outPre( argv[optind+1] );
    string outFQ1File = outPre + ".1.fq.gz";
    string outFQ2File = outPre + ".2.fq.gz";

    SeqLib::BamReader inBam;
    ifstream  inRdf(inRdfFile.c_str());
    ogzstream outFQ1(outFQ1File.c_str());
    ogzstream outFQ2(outFQ2File.c_str());

    if ( !inBam.Open(inBamFile) ) cerr << "open error: " << inBamFile  << endl, exit(1);
    if ( !inRdf.is_open()       ) cerr << "open error: " << inRdfFile  << endl, exit(1);

    if ( !outFQ1.rdbuf()->is_open() ) cerr << "open error: " << outFQ1File << endl, exit(1);
    if ( !outFQ2.rdbuf()->is_open() ) cerr << "open error: " << outFQ2File << endl, exit(1);

    SeqLib::BamWriter writer;
    if ( opt.outBamFile.size() > 0 ) {
        writer.Open(opt.outBamFile);
        if ( ! writer.IsOpen() ) cerr << "open error: " << opt.outBamFile << endl, exit(1);

        writer.SetHeader( removeReadGroup( inBam.Header() ) );
        writer.WriteHeader();
    }

    string famLine;
    size_t totalReadPair(0);

    if ( opt.randNread != 0 ) {
        while ( getline(inRdf, famLine) ) {
            istringstream stm(famLine);
            string chr, beg, end;
            size_t cnt;
            stm >> chr >> beg >> end >> cnt;
            totalReadPair += cnt/2;
        }
        inRdf.close();
        inRdf.open(inRdfFile.c_str());
    }

    if ( opt.randNread > totalReadPair )
        cerr << "WARNING: -n random select N read pairs > total " << totalReadPair << " read pairs in file";

    while ( getline(inRdf, famLine) ) {
        istringstream stm(famLine);
        string chr, beg, end;
        size_t cnt;
        stm >> chr >> beg >> end >> cnt;
        string chrBegEnd = chr + ":" + beg + "-" + end;

        mStrBrV watsonFam, crickFam;
        set<string> cigs;

        // read pairs
        for ( size_t i(0); i != cnt/2; i++ ) {
            BamRecord ra, rb;
            inBam.GetNextRecord(ra);
            inBam.GetNextRecord(rb);

            if ( ra.Qname() != rb.Qname() )
                cerr << "inBam and inBam.rdf is not consistent" << endl, exit(1);

            if ( ra.ReverseFlag() || !rb.ReverseFlag() )
                cerr << "not proper +/- pair: " << famLine << "\t" << i << endl, exit(1);

            if ( ra.MapQuality() < opt.mapQuaCutoff || rb.MapQuality() < opt.mapQuaCutoff )
                continue;

            string cigarAB = ra.CigarString() + ":" + rb.CigarString();

            // find S in cigar string
            if ( opt.filtSoftClip && cigarAB.find('S') != string::npos )
                continue;

            double rand_ratio = rand() / (double)RAND_MAX;
            if ( opt.randNread != 0
                    && rand_ratio > opt.randNread / (double)totalReadPair ) continue;

            if ( ra.FirstFlag() ) {
                watsonFam[cigarAB].push_back(ra);
                watsonFam[cigarAB].push_back(rb);
            }
            else {
                crickFam[cigarAB].push_back(ra);
                crickFam[cigarAB].push_back(rb);
            }

            cigs.insert(cigarAB);
        }

        // calibrate read family because cigar string might not be accurate
        calibrateFam(watsonFam);
        calibrateFam(crickFam);

        for ( auto &cg : cigs ) {
            mStrBrV::const_iterator w = watsonFam.find( cg );
            mStrBrV::const_iterator c =  crickFam.find( cg );

            if ( w != watsonFam.end() ) {
                if ( c != crickFam.end() ) {
                    if ( testCigarFam(watsonFam, crickFam, cg, opt) ) { // DCS
                        printConsensusRead(outFQ1, outFQ2, watsonFam, crickFam, cg, opt, chrBegEnd, writer );
                    }
                    else if (  opt.sscsOut
                            && w->second.size() <= opt.maxSupOnEachStrand * 2
                            && c->second.size() <= opt.maxSupOnEachStrand * 2 ) {
                        sscs(watsonFam, crickFam, cg, opt, chrBegEnd, writer);
                    }
                    else {
                        continue; // Fam is too big
                    }
                }
                else if ( opt.sscsOut && w->second.size() <= opt.maxSupOnEachStrand * 2 ) {
                    sscs(watsonFam, crickFam, cg, opt, chrBegEnd, writer);
                }
                else {
                    continue; // Fam is too big
                }
            }
            else if ( opt.sscsOut && c != crickFam.end() ) {
                if ( c->second.size() <= opt.maxSupOnEachStrand * 2 ) {
                    sscs(watsonFam, crickFam, cg, opt, chrBegEnd, writer);
                }
            }
            else {
                continue;
            }
        }
    }

    inBam.Close();
    inRdf.close();
    outFQ1.close();
    outFQ2.close();

    if ( opt.outBamFile.size() > 0 ) writer.Close();

    return 0;
}

