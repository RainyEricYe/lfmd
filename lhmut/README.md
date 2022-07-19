# Introduction
lhmut is a likelihood-based mutation detector for SNV and small indel.

# Preparatioon
alglib-3.16.0 (http://www.alglib.net/translator/re/alglib-3.16.0.cpp.gpl.tgz)

# Install
    git clone https://github.com/RainyEricYe/lhmut.git

    cd lhmut
    export ALGLIB=/path_to/alglib-3.16.0
    export LDPATH=/path_to_ld_path
    make

# Usage
    Program: lhmut (likelihood-based mutation detector)
    Version: v0.2
    Authors: Rui YE
    Contact: yerui@connect.hku.hk

    Options: lhmut -i in.pileup -o outfile

        -i [s]     input pileup file
        -o [s]     output mutation file
        -q [i]     base quality cutoff [20]
        -s [i]     min support num on each strand [1]
        -f [f]     min fraction of alterative allele in a read family [1e-5]
        -e [f]     precision of allele frequency while calculate likelihood ratio [1e-8]
        -g [f]     gap between likelihood ratios of major and rest genotypes [2.0]
        -x [i]     Encoding offset for phred quality scores [33]
        -d         debug mode [false]
        -h         help
        -v         version
