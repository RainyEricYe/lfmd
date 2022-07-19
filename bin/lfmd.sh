#!/bin/bash
set -euo pipefail

read inbam pre ref hash <<< "$@"

if [ $# -lt 3 ]; then
    echo "Usage: $0 in.sortByReadID.bam out_prefix ref"
    exit 1
fi

outd=`dirname $pre`
if [ ! -d $outd ]; then mkdir -p $outd; fi

bin=`dirname $0`
bamRdf_opt=""
bamDCS_opt=" -c -q 0 -b  "

echo `date` lfmd-start

# get read family
if [ ! -f $pre.family.bam.sign ]
then
    $bin/bamRdf $bamRdf_opt -i $inbam -o $pre.family.bam

    touch $pre.family.bam.sign
    echo `date` bamRdf-done
fi

# make SSCS or DCS
if [ ! -f $pre.dcs.bam.sign ]
then
    $bin/bamDCS $bamDCS_opt $pre.family.bam $pre.dcs -o $pre.dcs.bam

    touch $pre.dcs.bam.sign
    echo `date` bamDCS-done
fi

if [ ! -f $pre.dcs.sort.bam.sign ]
then
    samtools sort $pre.dcs.bam -o $pre.dcs.sort.bam
    samtools index $pre.dcs.sort.bam

    touch $pre.dcs.sort.bam.sign
    echo `date` sort-done
fi

if [ ! -f $pre.pileup.sign ]
then
    samtools mpileup -B -A -d 5000000 -f $ref $pre.dcs.sort.bam > $pre.pileup

    touch $pre.pileup.sign
    echo `date` pileup-done
fi

# call mutation
if [ ! -f $pre.mut.adj.sign ]
then
    $bin/lhmut -s 1 -f 0.00001 -e 1e-7 -i $pre.pileup -o $pre.mut
    $bin/adjustMut.pl $pre.mut | awk '$5+$10+$11>0' | grep -v F > $pre.mut.adj

    touch $pre.mut.adj.sign
    echo `date` lhmut-done
fi

echo `date` job-done
