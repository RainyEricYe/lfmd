# Introduction
bamRdf is a software which clusters reads of bam file into read families according to chromosome, start, and end position.

# Preparation
htslib(https://github.com/samtools/htslib)
SeqLib(https://github.com/walaj/SeqLib)

# Install
    git clone https://github.com/RainyEricYe/bamRdf.git

    cd bamRdf
    export SEQLIB=/path_to/SeqLib
    make

# Citation
    Rui Ye et al. LFMD: detecting low-frequency mutations in genome sequencing data without molecular tags.
    https://www.biorxiv.org/content/10.1101/617381v9

# Usage
    Program: bamRdf (cluster read pair by positon)
    Version: v1.4
    Authors: yerui
    Contact: yerui@connect.hku.hk

    Usage: bamRdf -i in.bam

       -t [s]    target region file (bed)
       -i [s]    input original bam file which sorted by Qname
       -o [s]    output bam file [out.bam]
       -d [i]    discard reads family which size < [int]  [0]
       -b [i]    memory control: N pairs of reads as a block [3000000]
       -m        mitochondrial mode. allow read start < 0 [false]
       -n [i]    limit number of read pairs [0, unlimited]
       -v        version
       -h        help

