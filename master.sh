#!/bin/bash
echo "Ok here we go"

source environ_qcd.sh

cat environ_qcd.sh

echo "Everything sourced"

# running pythia
python3 generate.py $1 $2 $3

hepmcfile=`ls *hepmc`
cp $hepmcfile ${ODIR_HEPMC}

# setting up delphes
git clone -b cl https://github.com/pumaphysics/delphes.git
cd delphes/
source /cvmfs/sft.cern.ch/lcg/views/LCG_92/x86_64-slc6-gcc62-opt/setup.sh; 
make -j1
sed -i 's/XXX/123/g' cards/papu_nopu/papu_CMS_PhaseII_HGCal.tcl

./DelphesHepMC3 cards/papu_nopu/papu_CMS_PhaseII_HGCal.tcl delphes.root $hepmcfile

./CLDelphes delphes.root flat.root

noext=""
noextfile="${hepmcfile/.hepmc/"$noext"}"  

#cp delphes.root ${ODIR_DELPHES}/${noextfile}_delphes.root
cp flat.root ${ODIR_FLAT}/${noextfile}_flat.root

