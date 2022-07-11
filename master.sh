#!/bin/bash
echo "Ok here we go"

source environ_qcd.sh

cat environ_qcd.sh

echo "Everything sourced"

python3 generate.py $1 $2 $3

cp *hepmc ${ODIR}

git clone -b papu https://github.com/pumaphysics/delphes.git

cd delphes/ 
source /cvmfs/sft.cern.ch/lcg/views/LCG_92/x86_64-slc6-gcc62-opt/setup.sh; 
make -j1



