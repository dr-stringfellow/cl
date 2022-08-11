import subprocess
import sys
import random
import uproot
import numpy as np
import awkward as ak
from os import path

print("#############################")
print("#############################")
print("######### Starting ##########")
print("#############################")
print("#############################")

# read branches
def read_branch(events,branchname):
    return events[branchname].array(entry_stop=1000000)

try:
    file = uproot.open(sys.argv[1])
except:
    exit(1)
events = file['events']
#print(events.keys())

minjetpt = 450

augmentations = open(sys.argv[2],'r').readlines() 
variations = []
for line in augmentations:
    v = line.strip().split(' ')[0]
    if v == 'nominal':
        continue
    variations.append(v)

print(variations)

parton_pt = read_branch(events,'parton_pt')
parton_eta = read_branch(events,'parton_eta')
parton_phi = read_branch(events,'parton_phi')
parton_e = read_branch(events,'parton_e')
jettype = read_branch(events,'jettype')
jet_pt = read_branch(events,'jet_pt')
jet_eta = read_branch(events,'jet_eta')
jet_phi = read_branch(events,'jet_phi')
jet_e = read_branch(events,'jet_e')
jet_msd = read_branch(events,'jet_msd')
jet_n2 = read_branch(events,'jet_n2')
pt = read_branch(events,'pt')
eta = read_branch(events,'eta')
phi = read_branch(events,'phi')
e = read_branch(events,'e')
charge = read_branch(events,'charge')
pdgid = read_branch(events,'pdgid')
d0 = read_branch(events,'d0')
dz = read_branch(events,'dz')

has_all_variations = []
variation_event_dict = {}

v_file = {}
v_events = {}
v_jettype = {}
v_parton_pt = {}
v_parton_eta = {}
v_parton_phi = {}
v_parton_e = {}
v_jet_pt = {}
v_jet_eta = {}
v_jet_phi = {}
v_jet_e = {}
v_jet_msd = {}
v_jet_n2 = {}
v_const_pt = {}
v_const_eta = {}
v_const_phi = {}
v_const_e = {}
v_const_charge = {}
v_const_pdgid = {}
v_const_d0 = {}
v_const_dz = {}

for v in variations:
    v_file[v] = uproot.open(sys.argv[1].replace("nominal",v))
    v_events[v] = v_file[v]['events']
    v_jettype[v] = read_branch(v_events[v],'jettype')
    v_parton_pt[v] = read_branch(v_events[v],'parton_pt')
    v_parton_eta[v] = read_branch(v_events[v],'parton_eta')
    v_parton_phi[v] = read_branch(v_events[v],'parton_phi')
    v_parton_e[v] = read_branch(v_events[v],'parton_e')
    v_jet_pt[v] = read_branch(v_events[v],'jet_pt')
    v_jet_eta[v] = read_branch(v_events[v],'jet_eta')
    v_jet_phi[v] = read_branch(v_events[v],'jet_phi')
    v_jet_e[v] = read_branch(v_events[v],'jet_e')
    v_jet_msd[v] = read_branch(v_events[v],'jet_msd')
    v_jet_n2[v] = read_branch(v_events[v],'jet_n2')
    v_const_pt[v] = read_branch(v_events[v],'pt')
    v_const_eta[v] = read_branch(v_events[v],'eta')
    v_const_phi[v] = read_branch(v_events[v],'phi')
    v_const_e[v] = read_branch(v_events[v],'e')
    v_const_charge[v] = read_branch(v_events[v],'charge')
    v_const_pdgid[v] = read_branch(v_events[v],'pdgid')
    v_const_d0[v] = read_branch(v_events[v],'d0')
    v_const_dz[v] = read_branch(v_events[v],'dz')

complete_counter = 0

dump_nominal = {}
dump_arrays = {}

for v in variations:
    dump_arrays[v] = {}

for evt in range(len(parton_pt)):

    if evt % 100 == 0:
        print(evt)

    #print(evt)
    #if evt > 30:
    #    break

    if jet_pt[evt] < minjetpt:
        continue
    all_checker = True

    #print("Nominal jet for event %i:"%evt, jet_pt[evt])

    #print(parton_pt[evt])
    #print(parton_eta[evt])
    #print(parton_phi[evt])
    #print(parton_e[evt])
    
    tmp_var_arrays = {}
    for v in variations:
        tmp_var_arrays[v] = {}

    for v in variations:

        mask = (np.array(v_parton_pt[v]) == parton_pt[evt]) & (np.array(v_parton_eta[v]) == parton_eta[evt]) & (np.array(v_parton_phi[v]) == parton_phi[evt]) & (np.array(v_parton_e[v]) == parton_e[evt])

        v_jettype_tmp = np.array(v_jettype[v])[mask]
        v_parton_pt_tmp = np.array(v_parton_pt[v])[mask]
        v_parton_eta_tmp = np.array(v_parton_eta[v])[mask]
        v_parton_phi_tmp = np.array(v_parton_phi[v])[mask]
        v_parton_e_tmp = np.array(v_parton_e[v])[mask]
        v_jet_pt_tmp = np.array(v_jet_pt[v])[mask]
        v_jet_eta_tmp = np.array(v_jet_eta[v])[mask]
        v_jet_phi_tmp = np.array(v_jet_phi[v])[mask]
        v_jet_e_tmp = np.array(v_jet_e[v])[mask]
        v_jet_msd_tmp = np.array(v_jet_msd[v])[mask]
        v_jet_n2_tmp = np.array(v_jet_n2[v])[mask]
        v_const_pt_tmp = np.array(v_const_pt[v])[mask]
        v_const_eta_tmp = np.array(v_const_eta[v])[mask]
        v_const_phi_tmp = np.array(v_const_phi[v])[mask]
        v_const_e_tmp = np.array(v_const_e[v])[mask]
        v_const_charge_tmp = np.array(v_const_charge[v])[mask]
        v_const_pdgid_tmp = np.array(v_const_pdgid[v])[mask]
        v_const_d0_tmp = np.array(v_const_d0[v])[mask]
        v_const_dz_tmp = np.array(v_const_dz[v])[mask]

        if len(v_parton_pt_tmp) != 1 or len(v_jet_pt_tmp) != 1 or v_jet_pt_tmp[0]<minjetpt:
            all_checker = False
            break            

        tmp_var_arrays[v]['jettype'] = v_jettype_tmp
        tmp_var_arrays[v]['parton_pt'] = v_parton_pt_tmp
        tmp_var_arrays[v]['parton_eta'] = v_parton_eta_tmp
        tmp_var_arrays[v]['parton_phi'] = v_parton_phi_tmp
        tmp_var_arrays[v]['parton_e'] = v_parton_e_tmp
        tmp_var_arrays[v]['jet_pt'] = v_jet_pt_tmp
        tmp_var_arrays[v]['jet_eta'] = v_jet_eta_tmp
        tmp_var_arrays[v]['jet_phi'] = v_jet_phi_tmp
        tmp_var_arrays[v]['jet_e'] = v_jet_e_tmp
        tmp_var_arrays[v]['jet_msd'] = v_jet_msd_tmp
        tmp_var_arrays[v]['jet_n2'] = v_jet_n2_tmp
        tmp_var_arrays[v]['pt'] = np.where(v_const_pt_tmp>0,np.log(v_const_pt_tmp),0)
        tmp_var_arrays[v]['relpt'] = np.where(v_const_pt_tmp>0,np.log(v_const_pt_tmp/v_jet_pt_tmp),0)
        tmp_var_arrays[v]['eta'] = v_const_eta_tmp-v_jet_eta_tmp
        tmp_var_arrays[v]['eta'] = np.where(v_const_pt_tmp>0,tmp_var_arrays[v]['eta'],0)
        tmp_var_arrays[v]['phi'] = v_const_phi_tmp-v_jet_phi_tmp
        tmp_var_arrays[v]['phi'] = np.where(v_const_pt_tmp>0,tmp_var_arrays[v]['phi'],0)
        tmp_var_arrays[v]['phi'] = np.where((tmp_var_arrays[v]['phi']<np.pi),tmp_var_arrays[v]['phi'],tmp_var_arrays[v]['phi']-2*np.pi)
        tmp_var_arrays[v]['phi'] = np.where((tmp_var_arrays[v]['phi']>-1*np.pi),tmp_var_arrays[v]['phi'],tmp_var_arrays[v]['phi']+2*np.pi)
        tmp_var_arrays[v]['dr'] = np.sqrt(tmp_var_arrays[v]['eta']*tmp_var_arrays[v]['eta']+tmp_var_arrays[v]['phi']*tmp_var_arrays[v]['phi'])
        tmp_var_arrays[v]['e'] = np.where(v_const_pt_tmp>0,np.log(v_const_e_tmp),0)
        tmp_var_arrays[v]['rele'] = np.where(v_const_pt_tmp>0,np.log(v_const_e_tmp/v_jet_e_tmp),0)
        tmp_var_arrays[v]['charge'] = v_const_charge_tmp
        tmp_var_arrays[v]['pdgid'] = v_const_pdgid_tmp
        tmp_var_arrays[v]['d0'] = np.tanh(v_const_d0_tmp)
        tmp_var_arrays[v]['dz'] = np.tanh(v_const_dz_tmp)


        if v == variations[-1]:
            complete_counter += 1

    if not all_checker:
        continue

    if all_checker:
        
        if 'jettype' not in dump_nominal:
            dump_nominal['jettype'] = []
            dump_nominal['jettype'].append(jettype[evt])
            dump_nominal['parton_pt'] = []
            dump_nominal['parton_pt'].append(parton_pt[evt])
            dump_nominal['parton_eta'] = []
            dump_nominal['parton_eta'].append(parton_eta[evt])
            dump_nominal['parton_phi'] = []
            dump_nominal['parton_phi'].append(parton_phi[evt])
            dump_nominal['parton_e'] = []
            dump_nominal['parton_e'].append(parton_e[evt])
            dump_nominal['jet_pt'] = []
            dump_nominal['jet_pt'].append(jet_pt[evt])
            dump_nominal['jet_eta'] = []
            dump_nominal['jet_eta'].append(jet_eta[evt])
            dump_nominal['jet_phi'] = []
            dump_nominal['jet_phi'].append(jet_phi[evt])
            dump_nominal['jet_e'] = []
            dump_nominal['jet_e'].append(jet_e[evt])
            dump_nominal['jet_msd'] = []
            dump_nominal['jet_msd'].append(jet_msd[evt])
            dump_nominal['jet_n2'] = []
            dump_nominal['jet_n2'].append(jet_n2[evt])
            dump_nominal['pt'] = []
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(pt[evt])),0)
            dump_nominal['pt'].append(tmppp)
            dump_nominal['relpt'] = []
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(pt[evt])/jet_pt[evt]),0)
            dump_nominal['relpt'].append(tmppp)
            dump_nominal['eta'] = []
            tmp_deta = np.array(eta[evt])-jet_eta[evt]
            tmp_deta = np.where(np.array(pt[evt])>0,tmp_deta,0)
            dump_nominal['eta'].append(tmp_deta)
            dump_nominal['phi'] = []
            tmp_dphi = np.array(phi[evt])-jet_phi[evt]
            tmp_dphi = np.where((tmp_dphi<np.pi),tmp_dphi,tmp_dphi-2*np.pi)
            tmp_dphi = np.where((tmp_dphi>-np.pi),tmp_dphi,tmp_dphi+2*np.pi)
            tmp_dphi = np.where(np.array(pt[evt])>0,tmp_dphi,0)
            dump_nominal['phi'].append(tmp_dphi)
            dump_nominal['e'] = []
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(e[evt])),0)
            dump_nominal['e'].append(tmppp)
            dump_nominal['rele'] = []
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(e[evt])/jet_e[evt]),0)
            dump_nominal['rele'].append(tmppp)
            dump_nominal['charge'] = []
            dump_nominal['charge'].append(np.array(charge[evt]))
            dump_nominal['pdgid'] = []
            dump_nominal['pdgid'].append(np.array(pdgid[evt]))
            dump_nominal['d0'] = []
            dump_nominal['d0'].append(np.tanh(np.array(d0[evt])))
            dump_nominal['dz'] = []
            dump_nominal['dz'].append(np.tanh(np.array(dz[evt])))
            dump_nominal['dr'] = []
            tmp_dr = np.sqrt(tmp_dphi*tmp_dphi+tmp_deta*tmp_deta)
            tmp_dr = np.where(np.array(pt[evt])>0,tmp_dr,0)
            dump_nominal['dr'].append(tmp_dr)
        else:
            dump_nominal['jettype'].append(jettype[evt])
            dump_nominal['parton_pt'].append(parton_pt[evt])
            dump_nominal['parton_eta'].append(parton_eta[evt])
            dump_nominal['parton_phi'].append(parton_phi[evt])
            dump_nominal['parton_e'].append(parton_e[evt])
            dump_nominal['jet_pt'].append(jet_pt[evt])
            dump_nominal['jet_eta'].append(jet_eta[evt])
            dump_nominal['jet_phi'].append(jet_phi[evt])
            dump_nominal['jet_e'].append(jet_e[evt])
            dump_nominal['jet_msd'].append(jet_msd[evt])
            dump_nominal['jet_n2'].append(jet_n2[evt])
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(pt[evt])),0)
            dump_nominal['pt'].append(tmppp)
            tmppp = np.where(np.array(pt[evt])>0,np.log(np.array(pt[evt])/jet_pt[evt]),0)
            dump_nominal['relpt'].append(tmppp)
            tmp_deta = np.array(eta[evt])-jet_eta[evt]
            tmp_deta = np.where(np.array(e[evt])>0,tmp_deta,0)
            dump_nominal['eta'].append(tmp_deta)
            tmp_dphi = np.array(phi[evt])-jet_phi[evt]
            tmp_dphi = np.where((tmp_dphi<np.pi),tmp_dphi,tmp_dphi-2*np.pi)
            tmp_dphi = np.where((tmp_dphi>-np.pi),tmp_dphi,tmp_dphi+2*np.pi)
            tmp_dphi = np.where(np.array(e[evt])>0,tmp_dphi,0)
            dump_nominal['phi'].append(tmp_dphi)
            tmppp = np.where(np.array(e[evt])>0,np.log(np.array(e[evt])),0)
            dump_nominal['e'].append(tmppp)
            tmppp = np.where(np.array(e[evt])>0,np.log(np.array(e[evt])/jet_e[evt]),0)
            dump_nominal['rele'].append(tmppp)
            dump_nominal['charge'].append(np.array(charge[evt]))
            dump_nominal['pdgid'].append(np.array(pdgid[evt]))
            dump_nominal['d0'].append(np.tanh(np.array(d0[evt])))
            dump_nominal['dz'].append(np.tanh(np.array(dz[evt])))
            tmp_dr = np.sqrt(tmp_dphi*tmp_dphi+tmp_deta*tmp_deta)
            tmp_dr = np.where(np.array(pt[evt])>0,tmp_dr,0)
            dump_nominal['dr'].append(tmp_dr)

        for v in variations:
            if 'jettype' not in dump_arrays[v]:
                dump_arrays[v]['jettype'] = []
                dump_arrays[v]['jettype'].append(tmp_var_arrays[v]['jettype'][0])
                dump_arrays[v]['parton_pt'] = []
                dump_arrays[v]['parton_pt'].append(tmp_var_arrays[v]['parton_pt'][0])
                dump_arrays[v]['parton_eta'] = []
                dump_arrays[v]['parton_eta'].append(tmp_var_arrays[v]['parton_eta'][0])
                dump_arrays[v]['parton_phi'] = []
                dump_arrays[v]['parton_phi'].append(tmp_var_arrays[v]['parton_phi'][0])
                dump_arrays[v]['parton_e'] = []
                dump_arrays[v]['parton_e'].append(tmp_var_arrays[v]['parton_e'][0])
                dump_arrays[v]['jet_pt'] = []
                dump_arrays[v]['jet_pt'].append(tmp_var_arrays[v]['jet_pt'][0])
                dump_arrays[v]['jet_eta'] = []
                dump_arrays[v]['jet_eta'].append(tmp_var_arrays[v]['jet_eta'][0])
                dump_arrays[v]['jet_phi'] = []
                dump_arrays[v]['jet_phi'].append(tmp_var_arrays[v]['jet_phi'][0])
                dump_arrays[v]['jet_e'] = []
                dump_arrays[v]['jet_e'].append(tmp_var_arrays[v]['jet_e'][0])
                dump_arrays[v]['jet_msd'] = []
                dump_arrays[v]['jet_msd'].append(tmp_var_arrays[v]['jet_msd'][0])
                dump_arrays[v]['jet_n2'] = []
                dump_arrays[v]['jet_n2'].append(tmp_var_arrays[v]['jet_n2'][0])
                dump_arrays[v]['pt'] = []
                dump_arrays[v]['pt'].append(tmp_var_arrays[v]['pt'][0])
                dump_arrays[v]['relpt'] = []
                dump_arrays[v]['relpt'].append(tmp_var_arrays[v]['relpt'][0])
                dump_arrays[v]['eta'] = []
                dump_arrays[v]['eta'].append(tmp_var_arrays[v]['eta'][0])
                dump_arrays[v]['phi'] = []
                dump_arrays[v]['phi'].append(tmp_var_arrays[v]['phi'][0])
                dump_arrays[v]['dr'] = []
                dump_arrays[v]['dr'].append(tmp_var_arrays[v]['dr'][0])
                dump_arrays[v]['e'] = []
                dump_arrays[v]['e'].append(tmp_var_arrays[v]['e'][0])
                dump_arrays[v]['rele'] = []
                dump_arrays[v]['rele'].append(tmp_var_arrays[v]['rele'][0])
                dump_arrays[v]['charge'] = []
                dump_arrays[v]['charge'].append(tmp_var_arrays[v]['charge'][0])
                dump_arrays[v]['pdgid'] = []
                dump_arrays[v]['pdgid'].append(tmp_var_arrays[v]['pdgid'][0])
                dump_arrays[v]['d0'] = []
                dump_arrays[v]['d0'].append(tmp_var_arrays[v]['d0'][0])
                dump_arrays[v]['dz'] = []
                dump_arrays[v]['dz'].append(tmp_var_arrays[v]['dz'][0])
            else:
                dump_arrays[v]['jettype'].append(tmp_var_arrays[v]['jettype'][0])
                dump_arrays[v]['parton_pt'].append(tmp_var_arrays[v]['parton_pt'][0])
                dump_arrays[v]['parton_eta'].append(tmp_var_arrays[v]['parton_eta'][0])
                dump_arrays[v]['parton_phi'].append(tmp_var_arrays[v]['parton_phi'][0])
                dump_arrays[v]['parton_e'].append(tmp_var_arrays[v]['parton_e'][0])
                dump_arrays[v]['jet_pt'].append(tmp_var_arrays[v]['jet_pt'][0])
                dump_arrays[v]['jet_eta'].append(tmp_var_arrays[v]['jet_eta'][0])
                dump_arrays[v]['jet_phi'].append(tmp_var_arrays[v]['jet_phi'][0])
                dump_arrays[v]['jet_e'].append(tmp_var_arrays[v]['jet_e'][0])
                dump_arrays[v]['jet_msd'].append(tmp_var_arrays[v]['jet_msd'][0])
                dump_arrays[v]['jet_n2'].append(tmp_var_arrays[v]['jet_n2'][0])
                dump_arrays[v]['pt'].append(tmp_var_arrays[v]['pt'][0])
                dump_arrays[v]['relpt'].append(tmp_var_arrays[v]['relpt'][0])
                dump_arrays[v]['eta'].append(tmp_var_arrays[v]['eta'][0])
                dump_arrays[v]['phi'].append(tmp_var_arrays[v]['phi'][0])
                dump_arrays[v]['dr'].append(tmp_var_arrays[v]['dr'][0])
                dump_arrays[v]['e'].append(tmp_var_arrays[v]['e'][0])
                dump_arrays[v]['rele'].append(tmp_var_arrays[v]['rele'][0])
                dump_arrays[v]['charge'].append(tmp_var_arrays[v]['charge'][0])
                dump_arrays[v]['pdgid'].append(tmp_var_arrays[v]['pdgid'][0])
                dump_arrays[v]['d0'].append(tmp_var_arrays[v]['d0'][0])
                dump_arrays[v]['dz'].append(tmp_var_arrays[v]['dz'][0])


rand_parton_pt = []
rand_parton_eta = []
rand_parton_phi = []
rand_parton_e = []
rand_jet_pt = []
rand_jet_eta = []
rand_jet_phi = []
rand_jet_e = []
rand_jet_msd = []
rand_jet_n2 = []
rand_pt = []
rand_relpt = []
rand_eta = []
rand_phi = []
rand_e = []
rand_rele = []
rand_charge = []
rand_pdgid = []
rand_d0 = []
rand_dz = []
rand_dr = []
rand_jettype = []
rand_vartype = []

for j in range(len(dump_nominal['jet_pt'])):
    idxrv = random.randrange(len(variations))
    rv = variations[idxrv]
    print("random_variation:", rv)
    rand_pt.append(dump_arrays[rv]['pt'][j])
    rand_relpt.append(dump_arrays[rv]['relpt'][j])
    rand_eta.append(dump_arrays[rv]['eta'][j])
    rand_phi.append(dump_arrays[rv]['phi'][j])
    rand_e.append(dump_arrays[rv]['e'][j])
    rand_rele.append(dump_arrays[rv]['rele'][j])
    rand_charge.append(dump_arrays[rv]['charge'][j])
    rand_pdgid.append(dump_arrays[rv]['pdgid'][j])
    rand_d0.append(dump_arrays[rv]['d0'][j])
    rand_dz.append(dump_arrays[rv]['dz'][j])
    rand_dr.append(dump_arrays[rv]['dr'][j])
    rand_parton_pt.append(dump_arrays[rv]['parton_pt'][j])
    rand_parton_eta.append(dump_arrays[rv]['parton_eta'][j])
    rand_parton_phi.append(dump_arrays[rv]['parton_phi'][j])
    rand_parton_e.append(dump_arrays[rv]['parton_e'][j])
    rand_jet_pt.append(dump_arrays[rv]['jet_pt'][j])
    rand_jet_eta.append(dump_arrays[rv]['jet_eta'][j])
    rand_jet_phi.append(dump_arrays[rv]['jet_phi'][j])
    rand_jet_e.append(dump_arrays[rv]['jet_e'][j])
    rand_jet_msd.append(dump_arrays[rv]['jet_msd'][j])
    rand_jet_n2.append(dump_arrays[rv]['jet_n2'][j])
    rand_jettype.append(dump_arrays[rv]['jettype'][j])
    rand_vartype.append(idxrv)


assert (np.array(rand_parton_eta) == np.array(dump_nominal['parton_eta'])).all()
for v in variations:
    assert (np.array(dump_arrays[v]['parton_eta']) == np.array(dump_nominal['parton_eta'])).all()


subprocess.call("rm -f %s"%sys.argv[1].replace('/flat','/skim').replace('nominal_flat.root',"random.root"),shell=True)
with uproot.recreate(sys.argv[1].replace('/flat','/skim').replace('nominal_flat.root',"random.root")) as ofile:
    ofile['events'] = {"pfcand": ak.zip({"pt": rand_pt, "relpt": rand_relpt, "eta": rand_eta, "phi": rand_phi, "e": rand_e, "charge": rand_charge, "pdgid": rand_pdgid, "d0": rand_d0, "dz": rand_dz, "dr": rand_dr, "rele": rand_rele}),"parton":ak.zip({"pt":rand_parton_pt,"eta":rand_parton_eta,"phi":rand_parton_phi,"e":rand_parton_e}),"jet":ak.zip({"pt":rand_jet_pt,"eta":rand_jet_eta,"phi":rand_jet_phi,"e":rand_jet_e,"msd":rand_jet_msd,"n2":rand_jet_n2}),'jettype':rand_jettype,'vartype':np.array(rand_vartype)}



#has_all_variations = True

#print(complete_counter,5000)

#print(dump_nominal['pt'])

print(variations)

for v in variations:
    #if v == "fsrRenHi2":
        #print("Here")
        #print(dump_arrays[v]['parton_pt'])

    assert (np.array(rand_parton_eta) == np.array(dump_nominal['parton_eta'])).all()


    if path.exists(sys.argv[1].replace('/flat','/skim').replace('nominal_flat.root',"%s.root"%v)):
        continue
    subprocess.call("rm -f %s"%sys.argv[1].replace('/flat','/skim').replace('nominal_flat.root',"%s.root"%v),shell=True)
    with uproot.recreate(sys.argv[1].replace('/flat','/skim').replace('nominal_flat.root',"%s.root"%v)) as ofile:
        ofile['events'] = {"pfcand": ak.zip({"pt": dump_arrays[v]['pt'], "relpt": dump_arrays[v]['relpt'], "eta": dump_arrays[v]['eta'], "phi": dump_arrays[v]['phi'], "e": dump_arrays[v]['e'], "charge": dump_arrays[v]['charge'], "pdgid": dump_arrays[v]['pdgid'], "d0": dump_arrays[v]['d0'], "dz": dump_arrays[v]['dz'], "dr": dump_arrays[v]['dr'], "rele": dump_arrays[v]['rele']}),"parton":ak.zip({"pt":dump_arrays[v]['parton_pt'],"eta":dump_arrays[v]['parton_eta'],"phi":dump_arrays[v]['parton_phi'],"e":dump_arrays[v]['parton_e']}),"jet":ak.zip({"pt":dump_arrays[v]['jet_pt'],"eta":dump_arrays[v]['jet_eta'],"phi":dump_arrays[v]['jet_phi'],"e":dump_arrays[v]['jet_e'],"msd":dump_arrays[v]['jet_msd'],"n2":dump_arrays[v]['jet_n2']}),'jettype':dump_arrays[v]['jettype']}

#if path.exists(sys.argv[1].replace('/flat','/skim').replace('_flat.root','.root')):
#    exit(1)

we = 4

#print(dump_nominal['eta'][we])
'''
print(dump_nominal['parton_pt'][we])
print(dump_nominal['parton_eta'][we])
print(dump_nominal['parton_phi'][we])
print(dump_nominal['jet_eta'][we])
print(dump_arrays['fsrRenLo']['parton_pt'][we])
print(dump_arrays['fsrRenLo']['parton_eta'][we])
print(dump_arrays['fsrRenLo']['parton_phi'][we])
#print(dump_arrays['fsrRenLo']['eta'][we])
print(dump_arrays['fsrRenLo']['jet_eta'][we])
'''

#print(dump_nominal['parton_pt'][2] == dump_arrays['fsrRenLo']['parton_pt'][2])

subprocess.call("rm -f %s"%sys.argv[1].replace('/flat','/skim').replace('_flat.root','.root'),shell=True)
with uproot.recreate(sys.argv[1].replace('/flat','/skim').replace('_flat.root','.root')) as ofile:
    ofile['events'] = {"pfcand": ak.zip({"pt": dump_nominal['pt'], "relpt": dump_nominal['relpt'], "eta": dump_nominal['eta'], "phi": dump_nominal['phi'], "e": dump_nominal['e'], "charge": dump_nominal['charge'], "pdgid": dump_nominal['pdgid'], "d0": dump_nominal['d0'], "dz": dump_nominal['dz'], "dr": dump_nominal['dr'], "rele": dump_nominal['rele']}), "parton": ak.zip({'pt':dump_nominal['parton_pt'],'eta':dump_nominal['parton_eta'],'phi':dump_nominal['parton_phi'],'e':dump_nominal['parton_e']}),"jet": ak.zip({'pt':dump_nominal['jet_pt'],'eta':dump_nominal['jet_eta'],'phi':dump_nominal['jet_phi'],'e':dump_nominal['jet_e'],'msd':dump_nominal['jet_msd'],'n2':dump_nominal['jet_n2']}),'jettype':dump_nominal['jettype']}
