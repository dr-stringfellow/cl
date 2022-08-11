#!/bin/bash
echo "Ok here we go"

environ=`ls environ*.sh`

source $environ

cat $environ

echo "Everything sourced"

# setting up Herwig7
source /afs/cern.ch/user/b/bmaier/herwig7/setup.sh

sed -i "s|XXXLHEFILEXXX|$1|g" lhe_template.in 

/cvmfs/pheno.egi.eu/Herwig/Herwig-7-2/bin/Herwig build lhe_template.in
/cvmfs/pheno.egi.eu/Herwig/Herwig-7-2/bin/Herwig run tmpprocess.run -N 20000

hepmcfile=`ls *hepmc`
#cp $hepmcfile ${ODIR_HEPMC}

# setting up delphes
git clone -b clhw7 https://github.com/pumaphysics/delphes.git
cd delphes/
source /cvmfs/sft.cern.ch/lcg/views/LCG_92/x86_64-slc6-gcc62-opt/setup.sh; 
make -j1
sed -i 's/XXX/123/g' cards/papu_nopu/papu_CMS_PhaseII_HGCal.tcl

./DelphesHepMC cards/papu_nopu/papu_CMS_PhaseII_HGCal.tcl delphes.root ../$hepmcfile

./CLDelphes delphes.root flat.root

infile=$1
base=`basename $infile`
noextfile="${base/.lhe/}"

#cp delphes.root ${ODIR_DELPHES}/${noextfile}_delphes.root
cp flat.root ${ODIR_FLAT}/${noextfile}_herwig_flat.root


