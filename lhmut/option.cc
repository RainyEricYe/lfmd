#include <iostream>
#include <string>
#include <unistd.h>

#include "option.h"
using namespace std;

void usage() {
    cout <<
        "Program: " << PROGRAM << " " << REMARKS << "\n"
        "Version: " << VERSION << "\n"
        "Authors: " << AUTHORS << "\n"
        "Contact: " << CONTACT << "\n\n"

        "Options: " << PROGRAM << " -i in.pileup -o outfile\n\n"
        "    -i [s]     input pileup file\n"
        "    -o [s]     output mutation file\n"

        "    -q [i]     base quality cutoff [20]\n"
        "    -s [i]     min support num on each strand [1]\n"
        "    -f [f]     min fraction of alterative allele in a read family [1e-5]\n"

        "    -e [f]     precision of allele frequency while calculate likelihood ratio [1e-8]\n"
        "    -g [f]     gap between likelihood ratios of major and rest genotypes [2.0]\n"
        "    -x [i]     Encoding offset for phred quality scores [33]\n"

        "    -d         debug mode [false]\n"
        "    -h         help\n"
        "    -v         version\n"
        "\n";
}

void Option::parse(int argc, char **argv)
{

    if ( argc < 2 ) usage(), exit(1);

    int c;
    while ( (c=getopt(argc,argv,"i:o:q:s:f:e:g:x:dvh")) != -1 ) {
        switch (c) {
            case 'i': infileName    = optarg;               break;
            case 'o': outfileName   = optarg;               break;

            case 'q': baseQuaCutoff      = atoi(optarg);    break;
            case 's': minSupOnEachStrand = atoi(optarg);    break;
            case 'f': minFractionInFam   = atof(optarg);    break;

            case 'e': freqPrecision = atof(optarg);         break;
            case 'g': lhrGapCutoff  = atof(optarg);         break;
            case 'x': phredOffset   = atoi(optarg);         break;

            case 'd': debug = true;                         break;
            case 'v': cerr << VERSION << endl;              exit(1);
            case 'h':
            default:  usage();                              exit(1);
        }
    }

    return;
}
