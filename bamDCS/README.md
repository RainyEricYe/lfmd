# Introduction
bamDCS is a software which generates double strand consensus reads for low-frequency mutation. The input files are from bamRdf.

# Preparation

* htslib(https://github.com/samtools/htslib)
- SeqLib(https://github.com/walaj/SeqLib)
* alglib(http://www.alglib.net/translator/re/alglib-3.16.0.cpp.gpl.tgz)

# Install
    git clone https://github.com/RainyEricYe/bamDCS.git

    cd bamDCS
    export SEQLIB=/path_to/SeqLib
    export ALGLIB=/path_to/alglib-3.16.0
    export LDPATH=/path_to_ld_path
    make

# Citation
    Rui Ye et al. LFMD: detecting low-frequency mutations in genome sequencing data without molecular tags.
    https://www.biorxiv.org/content/10.1101/617381v9
    
# Usage
    Program: bamDCS  (generate double strand consensus reads for low-frequency mutation)
    Version: v2.4
    Authors: yerui
    Contact: yerui@connect.hku.hk

    Options: bamDCS in.bam outprefix

        -q [i]     base quality cutoff [20]
        -Q [i]     map quality cutoff [30]
        -s [i]     min support num on each strand [3]
        -S [i]     max support num on each strand [3000]
        -N [f]     max fraction of N on consensus read [0.1]
        -f [f]     min fraction of alterative allele in a read family [0.002]
        -e [f]     precision of allele frequency while calculate likelihood ratio [0.00001]
        -g [f]     gap between likelihood ratios of major and rest genotypes [2.0]
        -x [i]     Encoding offset for phred quality scores [33]
        -t [i]     min support num to construct a haplotype seq [3]
        -c         discard soft-clipping reads [false]
        -C [i]     soft trim N base from both ends of read [5]
        -p [f]     expected PCR error rate during library construction [1e-5]
        -n [i]     N read pairs to be randomly selected and used [0, total]
        -o [s]     output bam File directly []
        -a         output single pair-end reads [false]
        -b         output SSCS [false]
        -d         debug mode [false]
        -h         help
        -v         version

