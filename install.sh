#!/bin/bash
set -euxo pipefail

cd `dirname $0`

# please modify the paths below
export SEQLIB=/home/yerui/work/src/SeqLib
export ALGLIB=/home/yerui/src/alglib-3.16.0
export LDPATH=/home/yerui/work/anaconda3

cd bamRdf; make; cd ..
cd bamDCS; make; cd ..
cd lhmut;  make; cd ..

cd bin
ln -sf ../bamRdf/bamRdf
ln -sf ../bamDCS/bamDCS
ln -sf ../lhmut/lhmut

echo job-done
