import numpy as np
import torch
from torch import nn
import h5py
from TorClasses import *
from contractions_handler import *
from datetime import datetime


# Save the results of Wicktrackt into an hdf5.

def flavor_specifier(flavor):
    if flavor in ['u', 'd', 'uB', 'dB', 'Y', 'YB']:
        return [7,7,7]#Light
    elif flavor in ['s', 'sB']:
        return [8,8,8]#Strange
    elif flavor in ['c', 'cB']:
        return [9,9,9]#Charm
    else:
        raise ValueError('Error in extracting quark flavor for wictrackt_to_np_extractor')

def flavorInsertion_specifier(flavor):
    if flavor in ['u', 'd', 'uB', 'dB', 'Y', 'YB']:
        return [17,17,17]#Light
    elif flavor in ['s', 'sB']:
        return [18,18,18]#Strange
    elif flavor in ['c', 'cB']:
        return [19,19,19]#Charm
    else:
        raise ValueError('Error in extracting quark flavor for wictrackt_to_np_extractor')

def Perambulator_Flavor(flavor):
    if np.all(flavor == np.array([7,7,7])):
        return 'Light'
    elif np.all(flavor == np.array([8,8,8])):
        return 'Strange'
    elif np.all(flavor == np.array([9,9,9])):
        return 'Charm'
    elif np.all(flavor == np.array([17,17,17])):
        return 'Light_Ins'
    elif np.all(flavor == np.array([18,18,18])):
        return 'Strange_Ins'
    elif np.all(flavor == np.array([19,19,19])):
        return 'Charm_Ins'
    else:
        raise ValueError('Error in extracting quark flavor for wictrackt_to_np_extractor') 

def wictrackt_to_np_extractor(dgrms_list, SumJ = False):
    list_of_diagrams = []
    hadron_seen      = []
    quark_seen       = set()
    for diagram in dgrms_list:
        if not SumJ:
            list_of_propagators = []
            for propagator in diagram.gpropagators():
                hadron1, hadron2 = propagator.ghdrn1(), propagator.ghdrn2()
                if hadron1 not in hadron_seen:
                    hadron_seen.append(hadron1)
                if hadron2 not in hadron_seen:
                    hadron_seen.append(hadron2)
                q1     = propagator.gbar()
                q2     = propagator.gnbar()
                flavor = flavor_specifier(propagator.gnbar().gflvr())
                Q_NB   = [q2.gtm(), q2.ghdrn_n(), q2.gqrk_hdrn_p()]
                Q_B    = [q1.gtm(), q1.ghdrn_n(), q1.gqrk_hdrn_p()]
                list_of_propagators.append([Q_NB, Q_B, flavor])
            list_of_diagrams.append([list_of_propagators, diagram.gff()])
        else:
            skip_diagram = False
            list_of_propagators = []
            list_of_propagators_conn_to_J = []
            for propagator in diagram.gpropagators():
                hadron1, hadron2 = propagator.ghdrn1(), propagator.ghdrn2()
                if hadron1 not in hadron_seen:
                    hadron_seen.append(hadron1)
                if hadron2 not in hadron_seen:
                    hadron_seen.append(hadron2)
                q1     = propagator.gbar()
                q2     = propagator.gnbar()

                # if the propagator is not connected to a current: business as usual
                #                   is connected to to a current only at one end: push propagator to current summation
                #                   is connected to a current at both ends: diagrams are not promoted to list_of_diagrams 
                if (q2.gtm() in [2,3]) or  (q1.gtm() in [2,3]):
                    if (q2.gtm()==2 and q1.gtm()==3) or (q2.gtm()==3 and q1.gtm()==2) or (q2.gtm()==q1.gtm()):
                        skip_diagram = True
                        break
                    else: 
                        list_of_propagators_conn_to_J.append(propagator)
                
                elif skip_diagram:
                    continue

                else:
                    flavor = flavor_specifier(propagator.gnbar().gflvr())
                    Q_NB   = [q2.gtm(), q2.ghdrn_n(), q2.gqrk_hdrn_p()]
                    Q_B    = [q1.gtm(), q1.ghdrn_n(), q1.gqrk_hdrn_p()]
                    print([Q_NB, Q_B, flavor]," gets added to propagators the normal way.")
                    list_of_propagators.append([Q_NB, Q_B, flavor])
                
                
            if (len(list_of_propagators_conn_to_J)!=0):
                if (len(list_of_propagators_conn_to_J) % 2 ==1):
                    print('Error: Odd number of propagators connected to currents')
                    continue
                # loop through all possible propagator pairs
                # only propagators with the same flavour can connect
                # remove redundant diagrams ( emerge when the same sink and source pairs can be connected via multiple currents
                for i in range(len(list_of_propagators_conn_to_J)):
                    for j in range(i+1, len(list_of_propagators_conn_to_J)):

                        
                        p1 = list_of_propagators_conn_to_J[i]
                        p2 = list_of_propagators_conn_to_J[j]
                        
                        q11 = p1.gbar()
                        q12 = p1.gnbar()
                        q21 = p2.gbar()
                        q22 = p2.gnbar()

                        if flavor_specifier(q11.gflvr())!=flavor_specifier(q21.gflvr()):
                            continue 
                
                        if q11.gtm() in [2, 3]:
                            flavor = flavorInsertion_specifier(q12.gflvr())
                            Q_NB   = [q12.gtm(), q12.ghdrn_n(), q12.gqrk_hdrn_p()]
                            Q_B    = [q21.gtm(), q21.ghdrn_n(), q21.gqrk_hdrn_p()]
                        else:
                            flavor = flavorInsertion_specifier(q12.gflvr())
                            Q_NB   = [q22.gtm(), q22.ghdrn_n(), q22.gqrk_hdrn_p()]
                            Q_B    = [q11.gtm(), q11.ghdrn_n(), q11.gqrk_hdrn_p()]

                        if [Q_NB, Q_B, flavor] not in list_of_propagators:
                            list_of_propagators.append([Q_NB, Q_B, flavor])
            if not any(prop_lst[0] == list_of_propagators for prop_lst in list_of_diagrams):  
                list_of_diagrams.append([list_of_propagators, diagram.gff()])
            #probably not necessary since the currents for mesons consist of 2 grassman variables so switching the currents shouldnt change gff
            #for this it is necessary to seperate the same flavour onto different currents
    if len(list_of_diagrams)!=0:
    	for i in range(0,len(list_of_diagrams)):
    		print('Afterwards the diagrams are:', list_of_diagrams[i][0])
    return list_of_diagrams


def writeresults(dgrms_list, *operators, SumJ = False, path = None):
    ymdhms = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path_T = f"contractions_list_{ymdhms}.hdf5"
    if (path is not None) and (not isinstance(path, str)):
        operators = (path,) + operators
        path = path_T
    if path is not None and (path.split('.')[-1] not in ['hdf5', 'h5']):
        print('Error in File name')
        print(f'File name is changed to: {path_T}')
        path = path_T
    if any(not isinstance(operator, OpTimeSlice) for operator in operators):
        raise TypeError('Error 0: For the writer method, operators must be defined on a time slice.')
    saving_layer = True
    counter = {0: 0, 1: 0, 2: 0, 3: 0}
    time_slice = {0: {}, 1: {}, 2: {}, 3: {}}
    for Operator in operators:
        O = Operator.hdrn_state
        t = Operator.tm
        if t not in [0, 1, 2, 3]:
            print("For the writer method, operators must be defined on a time slice.")
            print("The time must be one either: 0 (src) or 1 (snk) or 2 (current1) or 3 (current3).")
            print("The results of the diagram (in list form) will be saved without any information about the participating hadrons.")
            saving_layer = False
            break
        counter[t] += 1
        if isinstance(O, hdrn):
            time_slice[t] = O.ghdrn_t() + '('+ ''.join(q for q in O.gqrks()) + ')'
        elif isinstance(O, mhdrn):
            time_slice[t] = '[ '
            A_O = O.gAhdrn()
            for i, H in enumerate(A_O):
                if (i+1) != len(A_O):
                    p_k = '+ '
                else:
                    p_k = ' ]'
                time_slice[t] += H.ghdrn_t() + '('+ ''.join(q for q in H.gqrks()) + ') ' + p_k
        elif isinstance(O, (h2, h3, h4)):
            A_H = O.gAhdrn()
            time_slice[t] = ''.join(H.ghdrn_t() + '('+ ''.join(q for q in H.gqrks()) + ') ' for H in A_H)
        elif isinstance(O, (mh2, mh3, mh4)):
            time_slice[t] = '[ '
            A_O = O.gAhdrn()
            for i, A_H0 in enumerate(A_O):
                A_H = A_H0.gAhdrn()
                if (i+1) != len(A_O):
                    p_k = '+ '
                else:
                    p_k = ' ]'
                time_slice[t] += ''.join(H.ghdrn_t() + '('+ ''.join(q for q in H.gqrks()) + ') ' for H in A_H) + p_k
    list_of_diagrams = wictrackt_to_np_extractor(dgrms_list, SumJ = SumJ)
    with h5py.File(path, 'w') as write_data:
        if saving_layer:
            hadron_information = write_data.create_group('Hadron_Information')
            time_slice_info = []
            for t in [0, 1, 2, 3]:
                if isinstance(time_slice[t], str) and time_slice[t]:
                    info_line = f"Time slice {t}: {time_slice[t]}"
                else:
                    info_line = f"Time slice {t}: No Hadrons"
                time_slice_info.append(info_line)
            time_slice_dataset = hadron_information.create_dataset(
                "time_slice_info",
                data = np.array(time_slice_info, dtype = h5py.special_dtype(vlen=str))
            )
        results = write_data.create_group('Results')
        results.attrs['number_of_diagrams'] = len(list_of_diagrams)
        for i, diagram in enumerate(list_of_diagrams):
            diagram_group = results.create_group(f'Diagram_{i}')
            topology = diagram[0]
            factor   = diagram[1]
            topology_array = np.array(topology)
            diagram_group.create_dataset('Topology', data=topology_array)#[:]->file['Results']['Diagram_0']['Topology'][:]
            diagram_group.create_dataset('Factor', data=factor)#([])->file['Results']['Diagram_0']['Factor'][()]
    print('Results saved successfully')



# rewrite the results of a process into unique clusters
# The final result of this is at least as good as the result of Wicktrackt.. but here it is more precise and analyzed! 
def cluster_extractor(Path_Diagrams= None, Sub_Path = None):
# p is of the form     array([[1, 0, 0], [1, 1, 2], [7,7,7]]), i.e. [Q_NB, Q_B, flavor]
# d is of the form array([[[1, 0, 0],[1, 1, 2]],
# [[0, 1, 2], [0, 0, 0]], ...])
    def perambulator_extractor(p):
        Q, Q_Bar = p[0], p[1]
        H, H_Bar = Q[:2], Q_Bar[:2]
        Flavor   = Perambulator_Flavor(p[2])
        return Perambulator(Q, Q_Bar, H, H_Bar, Flavor)
    def diagram_extractor(dgrmn, d):
        perambulators = []
        for p in d:
            perambulators.append(perambulator_extractor(p))
        return Diagram(dgrmn, *perambulators).organize()
    organized_diagrams            = []
    numerical_factors_of_diagrams = []
    with h5py.File(Path_Diagrams, 'r') as all_results0:
        if Sub_Path is not None:
            all_results = all_results0[Sub_Path]
        else:
            all_results = all_results0
        N = all_results['Results'].attrs['number_of_diagrams']
        for i in range(N):
            numerical_factors_of_diagrams.append(all_results['Results']['Diagram_'+str(i)]['Factor'][()])
            organized_diagrams.append(diagram_extractor(i, all_results['Results']['Diagram_'+str(i)]['Topology'][:]))
    all_results_unclustered = Diagram_Container(organized_diagrams)
    all_results_clusterd   = all_results_unclustered.do_clustering()
    print(f'Diagrams have been successfully clustered! There are {len(all_results_clusterd)} clusters')
    return all_results_clusterd, numerical_factors_of_diagrams
