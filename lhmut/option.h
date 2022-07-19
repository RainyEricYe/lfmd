#ifndef OPTION_H_
#define OPTION_H_

#define PROGRAM "lhmut"
#define VERSION "v0.2"
#define AUTHORS "Rui YE"
#define CONTACT "yerui@connect.hku.hk"
#define REMARKS "(likelihood-based mutation detector)"

#include <string>
using namespace std;

void usage();

class Option {
    public:
        Option():
            infileName(""),
            outfileName(""),
            baseQuaCutoff(20),
            minSupOnEachStrand(1),
            minFractionInFam(1e-5),
            freqPrecision(1e-8),
            lhrGapCutoff(2.0),
            phredOffset(33),
            debug(false) {}

        ~Option() {}

        void parse(int argc, char **argv);

        string infileName;
        string outfileName;
        int    baseQuaCutoff;
        unsigned long  minSupOnEachStrand;
        double minFractionInFam;
        double freqPrecision;
        double lhrGapCutoff;
        int    phredOffset;
        bool   debug;

};

#endif // OPTION_H_
