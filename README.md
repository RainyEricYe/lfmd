# LFMD
Low Frequency Mutation Detector

# Preparation

* htslib(https://github.com/samtools/htslib)
- SeqLib(https://github.com/walaj/SeqLib)
* alglib(http://www.alglib.net/translator/re/alglib-3.15.0.cpp.gpl.tgz)
- bamRdf(https://github.com/RainyEricYe/bamRdf)
* bamDCS(https://github.com/RainyEricYe/bamDCS)
- lhmut(https://github.com/RainyEricYe/lhmut)

# Install
    git clone https://github.com/RainyEricYe/lfmd2.git
    cd lfmd2

    cd example; sh test.sh # check if it works

    # if not work, you need to compile it as follows.
    cd ..
    vi install.sh
    # modify paths to SeqLib alglib and librarys    
    sh install.sh
    cd example; sh test.sh

    

# Contact
  yerui@connect.hku.hk
  
# Citation
Rui Ye et al. LFMD: detecting low-frequency mutations in genome sequencing data without molecular tags. https://www.biorxiv.org/content/10.1101/617381v11
  
# Usage
    vi lfmd.sh && ajust parameters and file paths
    ./lfmd.sh in.sortByReadID.bam outprefix
