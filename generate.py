from numpythia import Pythia, hepmc_write, hepmc_read
from numpythia import STATUS, HAS_END_VERTEX, ABS_PDG_ID
from numpythia.testcmnd import get_cmnd
from numpy.testing import assert_array_equal
import sys

lhefile = sys.argv[1]

update = True
for a in sys.argv:
    if 'json' in a:
        update = False
    
nickname = sys.argv[2]
setting = sys.argv[3]

key = setting.split('=')[0]
value = setting.split('=')[1]

print(nickname)
print(key)
print(value)


#pythia = Pythia(config=sys.argv[1], random_state=1)

#https://github.com/cms-sw/cmssw/blob/master/Configuration/Generator/python/Pythia8CommonSettings_cfi.py
common_dict = {
    'Tune:preferLHAPDF': 2,
    'Main:timesAllowErrors': 10000,
    'Check:epTolErr': 0.01,
    'Beams:setProductionScalesFromLHEF': 'off',
    'SLHA:minMassSM': 1000.,
    'ParticleDecays:limitTau0': 'on',
    'ParticleDecays:tau0Max': 10,
    'ParticleDecays:allowPhotonRadiation': 'on',
    '23:onMode': 'off',
    '23:onIfAny': 13,
}

#https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/Configuration/Generator/python/MCTunes2017/PythiaCP5Settings_cfi.py
cp5_dict = {
    'Tune:pp': 14,
    'Tune:ee': 7,
    'MultipartonInteractions:ecmPow': 0.03344,
    'MultipartonInteractions:bProfile': 2,
    'MultipartonInteractions:pT0Ref': 1.41,
    'MultipartonInteractions:coreRadius': 0.7634,
    'MultipartonInteractions:coreFraction': 0.63,
    'ColourReconnection:range': 5.176,
    'SigmaTotal:zeroAXB': 'off',
    'SpaceShower:alphaSorder': 2,
    'SpaceShower:alphaSvalue': 0.118,
    'SigmaProcess:alphaSvalue': 0.118,
    'SigmaProcess:alphaSorder': 2,
    'MultipartonInteractions:alphaSvalue': 0.118,
    'MultipartonInteractions:alphaSorder': 2,
    'TimeShower:alphaSorder': 2,
    'TimeShower:alphaSvalue': 0.118,
    'SigmaTotal:mode': 0,
    'SigmaTotal:sigmaEl': 21.89,
    'SigmaTotal:sigmaTot': 100.309,
    'PDF:pSet':20,
}

cp5down_dict = {
        'MultipartonInteractions:pT0Ref':1.46,
        'MultipartonInteractions:coreRadius':0.6879,
        'MultipartonInteractions:coreFraction':0.7253,
        'ColourReconnection:range':4.691,
}

cp5up_dict = {
        'MultipartonInteractions:pT0Ref':1.407,
        'MultipartonInteractions:coreRadius':0.6671,
        'MultipartonInteractions:coreFraction':0.4281,
        'ColourReconnection:range':4.881,
}

seeddown_dict = {
        'Tune:preferLHAPDF': 2,
}

seedup_dict = {
        'Tune:preferLHAPDF': 2,
}

params_dict = {'Beams:eCM' : 13000.,'Beams:frameType' : 4,'Beams:LHEF':lhefile}
params_dict.update(common_dict)
params_dict.update(cp5_dict)

if update and nickname != "nominal":
    params_dict.update({key:float(value)})
elif nickname == 'cp5up':
    params_dict.update(cp5up_dict)
elif nickname == 'cp5down':
    params_dict.update(cp5down_dict)
elif nickname == 'seedup':
    params_dict.update(seedup_dict)
elif nickname == 'seeddown':
    params_dict.update(seeddown_dict)


import random
randomseed = random.randrange(1e8)

pythia = Pythia(params=params_dict,random_state=randomseed)

selection = ((STATUS == 1) & ~HAS_END_VERTEX &\
            (ABS_PDG_ID != 12) & (ABS_PDG_ID != 14) & (ABS_PDG_ID != 16))


ofile = lhefile.split("/")[-1].replace('.lhe','_%s.hepmc'%nickname)
# generate events while writing to ascii hepmc
for event in hepmc_write(ofile, pythia(events=10000)):
    array1 = event.all(selection)

