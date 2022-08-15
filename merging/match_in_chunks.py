import h5py
import awkward as ak
import numpy as np
import uproot as uproot
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys
import subprocess

import tqdm

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

dijet_path = sys.argv[1]
higgs_path = sys.argv[2]
maxevents = 3000000

dijet_events = 0
higgs_events = 0

dijet_files = glob.glob('%s/*nominal.root'%dijet_path)
higgs_files = glob.glob('%s/*nominal.root'%higgs_path)



#list_dijet_files = chunks(dijet_files,20)
#list_higgs_files = chunks(higgs_files,10)

branches = ['pfcand_pt','pfcand_relpt','pfcand_eta','pfcand_phi','pfcand_e','pfcand_rele','pfcand_charge','pfcand_d0','pfcand_dz','pfcand_dr']
branches_pid = ['pfcand_pdgid']

parton_branches = ['parton_pt','parton_eta','parton_phi','parton_e']
jet_branches = ['jet_pt','jet_eta','jet_phi','jet_e','jet_msd','jet_n2']



def do_merge(qcd_chunk,higgs_chunk,cidx):

    print("Chunk %i"%cidx)
    
    def create_arrays(files):
        dijet_counter = 0
        large_dijet_features = None
        large_dijet_aug_features = None
        large_parton_dijet_features = None
        large_parton_dijet_aug_features = None
        large_jet_dijet_features = None
        large_jet_dijet_aug_features = None
        large_vartype = None
        large_vartype_aug = None
        large_all_photon = None
        large_all_electron = None
        large_all_muon = None
        large_all_chs = None
        large_all_nhs = None
        large_all_aug_photon = None
        large_all_aug_electron = None
        large_all_aug_muon = None
        large_all_aug_chs = None
        large_all_aug_nhs = None
        
        for f in tqdm.tqdm(files):
            
            print(f)
            
            try:
                #print("A")
                dijet_events = uproot.open(f)['events']
                #print("B")
                dijet_nevents = dijet_events.num_entries
                dijet_aug_events = uproot.open(f.replace("nominal","random"))['events']
                dijet_aug_nevents = dijet_aug_events.num_entries
                #print("C")
                if dijet_nevents != dijet_aug_nevents:
                    continue
            except:
                continue


            feature_array = np.stack([np.array(dijet_events[b].array()) for b in branches],axis=-1)
            feature_pid_array = np.stack([np.array(dijet_events[b].array()) for b in branches_pid],axis=-1)
            is_photon = 1*(feature_pid_array == 22)
            is_electron = 1*(np.abs(feature_pid_array) == 11)
            is_muon = 1*(np.abs(feature_pid_array) == 13)
            is_ch = 1* ( (np.abs(feature_pid_array) == 211) | (np.abs(feature_pid_array) == 321) | (np.abs(feature_pid_array) == 2212)  | (np.abs(feature_pid_array) == 3112) | (np.abs(feature_pid_array) == 3222) | (np.abs(feature_pid_array) == 3312) | (np.abs(feature_pid_array) == 3334)  )
            is_nh = 1*(feature_pid_array == 0)

            feature_aug_array = np.stack([np.array(dijet_aug_events[b].array()) for b in branches],axis=-1)
            feature_aug_pid_array = np.stack([np.array(dijet_aug_events[b].array()) for b in branches_pid],axis=-1)
            is_aug_photon = 1*(feature_aug_pid_array == 22)
            is_aug_electron = 1*(np.abs(feature_aug_pid_array) == 11)
            is_aug_muon = 1*(np.abs(feature_aug_pid_array) == 13)
            is_aug_ch = 1* ( (np.abs(feature_aug_pid_array) == 211) | (np.abs(feature_aug_pid_array) == 321) | (np.abs(feature_aug_pid_array) == 2212)  | (np.abs(feature_aug_pid_array) == 3112) | (np.abs(feature_aug_pid_array) == 3222) | (np.abs(feature_aug_pid_array) == 3312) | (np.abs(feature_aug_pid_array) == 3334)  )
            is_aug_nh = 1*(feature_aug_pid_array == 0)
            

            parton_feature_array = np.stack([np.array(dijet_events[b].array()) for b in parton_branches],axis=-1)
            parton_feature_aug_array = np.stack([np.array(dijet_aug_events[b].array()) for b in parton_branches],axis=-1)
            jet_feature_array = np.stack([np.array(dijet_events[b].array()) for b in jet_branches],axis=-1)
            jet_feature_aug_array = np.stack([np.array(dijet_aug_events[b].array()) for b in jet_branches],axis=-1)
            vartype_array = np.array([-1 for i in range(len(dijet_events['jettype'].array()))])
            vartype_aug_array = np.array(dijet_aug_events['vartype'].array())
            
            dijet_counter += dijet_nevents
            
            assert feature_array.shape == feature_aug_array.shape
            
            #print(feature_array.shape)
            #print(feature_array)
            
            if large_dijet_features is None:
                large_dijet_features = feature_array
                large_dijet_aug_features = feature_aug_array
                large_parton_dijet_features = parton_feature_array
                large_parton_dijet_aug_features = parton_feature_aug_array
                large_jet_dijet_features = jet_feature_array
                large_jet_dijet_aug_features = jet_feature_aug_array
                large_vartype = vartype_array
                large_vartype_aug = vartype_aug_array
                large_all_photon = is_photon
                large_all_electron = is_electron
                large_all_muon = is_muon
                large_all_chs = is_ch
                large_all_nhs = is_nh
                large_all_aug_photon = is_aug_photon
                large_all_aug_electron = is_aug_electron
                large_all_aug_muon = is_aug_muon
                large_all_aug_chs = is_aug_ch
                large_all_aug_nhs = is_aug_nh
            else:
                large_dijet_features = np.concatenate((large_dijet_features,feature_array))
                large_dijet_aug_features = np.concatenate((large_dijet_aug_features,feature_aug_array))
                large_parton_dijet_features = np.concatenate((large_parton_dijet_features,parton_feature_array))
                large_parton_dijet_aug_features = np.concatenate((large_parton_dijet_aug_features,parton_feature_aug_array))
                large_jet_dijet_features = np.concatenate((large_jet_dijet_features,jet_feature_array))
                large_jet_dijet_aug_features = np.concatenate((large_jet_dijet_aug_features,jet_feature_aug_array))
                large_vartype = np.concatenate((large_vartype,vartype_array))
                large_vartype_aug = np.concatenate((large_vartype_aug,vartype_aug_array))
                large_all_photon = np.concatenate((large_all_photon,is_photon))
                large_all_electron = np.concatenate((large_all_electron,is_electron))
                large_all_muon = np.concatenate((large_all_muon,is_muon))
                large_all_chs = np.concatenate((large_all_chs,is_ch))
                large_all_nhs = np.concatenate((large_all_nhs,is_nh))
                large_all_aug_photon = np.concatenate((large_all_aug_photon,is_aug_photon))
                large_all_aug_electron = np.concatenate((large_all_aug_electron,is_aug_electron))
                large_all_aug_muon = np.concatenate((large_all_aug_muon,is_aug_muon))
                large_all_aug_chs = np.concatenate((large_all_aug_chs,is_aug_ch))
                large_all_aug_nhs = np.concatenate((large_all_aug_nhs,is_aug_nh))
            
                #print(large_dijet_aug_features.shape)
                
                
            #print(float(dijet_counter/nevents))


    
        #print(large_dijet_features)
        #print(large_all_nhs)
        ffinal_dijet_features = np.concatenate((large_dijet_features,large_all_photon,large_all_muon,large_all_electron,large_all_chs,large_all_nhs),axis=-1)
        ffinal_dijet_aug_features = np.concatenate((large_dijet_aug_features,large_all_aug_photon,large_all_aug_muon,large_all_aug_electron,large_all_aug_chs,large_all_aug_nhs),axis=-1)
        #print(final_dijet_features.shape)
        #print("A")
        #print(large_all_flat)
        #print("B")
        #print(large_all_flat.shape)
        #print("C")
        #for u in np.unique(large_all_flat):
        #    print(u)
        #print(np.unique(large_all_flat))
        
        
        #return large_dijet_features[:nevents], large_dijet_aug_features[:nevents], large_parton_dijet_features[:nevents], large_parton_dijet_aug_features[:nevents], large_jet_dijet_features[:nevents], large_jet_dijet_aug_features[:nevents], large_vartype, large_vartype_aug
        return ffinal_dijet_features, ffinal_dijet_aug_features, large_parton_dijet_features, large_parton_dijet_aug_features, large_jet_dijet_features, large_jet_dijet_aug_features, large_vartype, large_vartype_aug


    dijet_features, dijet_aug_features, parton_dijet_features, parton_dijet_aug_features, jet_dijet_features, jet_dijet_aug_features, dijet_vartype, dijet_vartype_aug = create_arrays(qcd_chunk)
    higgs_features, higgs_aug_features, parton_higgs_features, parton_higgs_aug_features, jet_higgs_features, jet_higgs_aug_features, higgs_vartype, higgs_vartype_aug = create_arrays(higgs_chunk)


    min_events = min(len(dijet_features),len(higgs_features))

    #print("Min_events in chunk %i"%cidx)
    #print(min_events)
    
    dijet_features = dijet_features[:min_events]
    dijet_aug_features = dijet_aug_features[:min_events]
    higgs_features = higgs_features[:min_events]
    higgs_aug_features = higgs_aug_features[:min_events]
    parton_dijet_features = parton_dijet_features[:min_events]
    parton_dijet_aug_features = parton_dijet_aug_features[:min_events]
    parton_higgs_features = parton_higgs_features[:min_events]
    parton_higgs_aug_features = parton_higgs_aug_features[:min_events]
    jet_dijet_features = jet_dijet_features[:min_events]
    jet_dijet_aug_features = jet_dijet_aug_features[:min_events]
    jet_higgs_features = jet_higgs_features[:min_events]
    jet_higgs_aug_features = jet_higgs_aug_features[:min_events]
    dijet_vartype = dijet_vartype[:min_events]
    dijet_vartype_aug = dijet_vartype_aug[:min_events]
    higgs_vartype = higgs_vartype[:min_events]
    higgs_vartype_aug = higgs_vartype_aug[:min_events]
    
    truth_label = np.array([0 for i in range(len(dijet_features))] + [1 for i in range(len(higgs_features))])
    
    indices = [i for i in range(len(truth_label))]
    
    import random
    random.shuffle(indices)
    
    #print(indices)
    #print(truth_label)


    final_features = np.array(ak.Array(np.concatenate((dijet_features,higgs_features)))[indices])
    final_aug_features = np.array(ak.Array(np.concatenate((dijet_aug_features,higgs_aug_features)))[indices])
    final_indices = np.array(ak.Array(np.array(truth_label))[indices])
    
    final_parton_features = np.array(ak.Array(np.concatenate((parton_dijet_features,parton_higgs_features)))[indices])
    final_parton_aug_features = np.array(ak.Array(np.concatenate((parton_dijet_aug_features,parton_higgs_aug_features)))[indices])
    
    final_jet_features = np.array(ak.Array(np.concatenate((jet_dijet_features,jet_higgs_features)))[indices])
    final_jet_aug_features = np.array(ak.Array(np.concatenate((jet_dijet_aug_features,jet_higgs_aug_features)))[indices])
    
    final_vartype = np.array(ak.Array(np.concatenate((dijet_vartype,higgs_vartype)))[indices])
    final_vartype_aug = np.array(ak.Array(np.concatenate((dijet_vartype_aug,higgs_vartype_aug)))[indices])
    
    c = np.zeros(shape=(2*final_features.shape[0], final_features.shape[1],final_features.shape[2]))
    #c = np.empty((final_features.size + final_aug_features.size,), dtype=final_features.dtype)
    c[0::2] = final_features
    c[1::2] = final_aug_features
    
    assert (final_parton_features == final_parton_aug_features).all()
    
    parton_c = np.zeros(shape=(2*final_parton_features.shape[0], final_parton_features.shape[1]))
    parton_c[0::2] = final_parton_features
    parton_c[1::2] = final_parton_aug_features
    
    jet_c = np.zeros(shape=(2*final_jet_features.shape[0], final_jet_features.shape[1]))
    jet_c[0::2] = final_jet_features
    jet_c[1::2] = final_jet_aug_features
    
    #print(final_vartype)

    vartype_c = np.zeros(shape=(2*final_vartype.shape[0]))
    vartype_c[0::2] = final_vartype
    vartype_c[1::2] = final_vartype_aug
    
    #print(parton_c)
    
    idx = np.repeat(final_indices, 2)
    #print(idx)
    #print(final_features[0])
    #print(final_aug_features[0])
    #print(c[0])
    #print(c[1])
    #print(dijet_features.shape)
    #print(higgs_features.shape)



    def closestNumber(n, m) :
        # Find the quotient
        q = int(n / m)
     
        # 1st possible closest number
        n1 = m * q
        
        # 2nd possible closest number
        if((n * m) > 0) :
            n2 = (m * (q + 1))
        else :
            n2 = (m * (q - 1))
            
        # if true, then n1 is the required closest number
        if (abs(n - n1) < abs(n - n2)) :
            return n1
            
        # else n2 is the required closest number
        return n2


    #print("How long is the parton_c array?")
    #print(len(parton_c))

    closest_divisor = closestNumber(len(parton_c),3)
    if closest_divisor % 2 != 0:
        closest_divisor -= 3

    one_third_events = closest_divisor/3
    
    for i,j in enumerate(["train","test","val"]):
        #print("A")
        #print(i*one_third_events)
        #print((i+1)*one_third_events)
        hf = h5py.File(sys.argv[3]+"/%s/"%j+"merged_%i_%s.h5"%(cidx,j), 'w')
        #hf.create_dataset('features', data=c[int(0.333*(i)*len(c)):int(0.333*(i+1)*len(c))])
        #print("Asserting for chunk %i"%cidx)
        #print(len(c[int(i*one_third_events):int((i+1)*one_third_events)]) % 2)
        assert len(c[int(i*one_third_events):int((i+1)*one_third_events)]) % 2 == 0
        assert (parton_c[int(i*one_third_events):int((i+1)*one_third_events)][0::2] == parton_c[int(i*one_third_events):int((i+1)*one_third_events)][1::2]).all()
        hf.create_dataset('features', data=c[int(i*one_third_events):int((i+1)*one_third_events)])
        #hf.create_dataset('parton_features', data=parton_c[int(0.333*(i)*len(parton_c)):int(0.333*(i+1)*len(parton_c))])
        hf.create_dataset('parton_features', data=parton_c[int(i*one_third_events):int((i+1)*one_third_events)])
        #hf.create_dataset('jet_features', data=jet_c[int(0.333*(i)*len(jet_c)):int(0.333*(i+1)*len(jet_c))])
        hf.create_dataset('jet_features', data=jet_c[int(i*one_third_events):int((i+1)*one_third_events)])
        #hf.create_dataset('truth_label', data=idx[int(0.333*(i)*len(idx)):int(0.333*(i+1)*len(idx))])
        hf.create_dataset('truth_label', data=idx[int(i*one_third_events):int((i+1)*one_third_events)])
        #hf.create_dataset('vartype', data=vartype_c[int(0.333*(i)*len(vartype_c)):int(0.333*(i+1)*len(vartype_c))])
        hf.create_dataset('vartype', data=vartype_c[int(i*one_third_events):int((i+1)*one_third_events)])
        hf.close()
    
        






higgs_file_idx = 0
qcd_files_in_chunks = []
higgs_files_in_chunks = []
total_qcd_events = 0
for i,c in enumerate(chunks(dijet_files,10)):
    print("Chunk ... %i"%i)
    tmp_qcdfiles = []
    events_in_chunk = 0    
    for cf in c:
        tmp_qcdfiles.append(cf)
        #print(cf)
        cftree = uproot.open(cf)['events']
        events_in_chunk += cftree.num_entries
    print("Number of events %i"%events_in_chunk)
    qcd_files_in_chunks.append(tmp_qcdfiles)
    total_qcd_events += events_in_chunk
    
    higgs_events_in_chunk = 0 
    higgs_files_in_chunk = []
    while higgs_file_idx < len(higgs_files) and higgs_events_in_chunk < events_in_chunk:
        #print(higgs_files[higgs_file_idx])
        higgstree = uproot.open(higgs_files[higgs_file_idx])['events']
        higgs_events_in_chunk += higgstree.num_entries
        higgs_files_in_chunk.append(higgs_files[higgs_file_idx])
        higgs_file_idx += 1

    higgs_files_in_chunks.append(higgs_files_in_chunk)

    if total_qcd_events > float(maxevents)/2.:
        break



for i in range(len(qcd_files_in_chunks)):
    if i == len(qcd_files_in_chunks):
        break
    do_merge(qcd_files_in_chunks[i],higgs_files_in_chunks[i],i)
    
exit(0)
