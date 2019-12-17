import sys
#add local paths to load custom methods
if not any(('scripts' in pth for pth in sys.path)):
    sys.path.append('scripts/')


from interface_methods import *
from interface_formation import formTime
from itertools import product
from collections import defaultdict
import pandas as pd
import random

def readEvolutionRecord(mu,S_c,rate,duplicate=True,FNAME='/rscratch/asl47/Duplication/EvoRecords/'):

    simulations, sets = [[]], []

    for line in open(f'{FNAME}EvoRecord_Mu{mu:.6f}_S{S_c:.6f}_D{rate:.6f}.txt'):
        #blank link, add new entry if previous entry has data
        if line == '\n':
            if simulations[-1]:
                simulations.append([])
        #information on phenotype evolution order
        elif ',' in line:
            sets.append(getSuperSets([{int(i) for i in ordering.split()} for ordering in line.rstrip().split(',') if ordering]))
          #information on phenotype polyomino
        else:
            parts=line.rstrip().split()
            simulations[-1].append(tuple(int(i) for i in parts[:4])+tuple(tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)))

    return simulations,sets


def getSuperSets(list_of_sets):
    super_sets=[]
    for l1 in range(len(list_of_sets)):
        for l2 in range(l1+1,len(list_of_sets)):
            if list_of_sets[l1] < list_of_sets[l2]:
                break
        else:
            super_sets.append(sorted(list_of_sets[l1]))
    return super_sets
 
def cleanRecord(full_simulations,full_sets):

    def timeCondition(phenotypes, sequence):
        if len(sequence) == 1:
            return sequence
        if not sequence or len(sequence) < 2:
            return
        new_sequence = sequence[:]

        for i in range(len(sequence)-1):
            if phenotypes[sequence[i]][2]==phenotypes[sequence[i+1]][2]:
                if phenotypes[sequence[i]][0]>=phenotypes[sequence[i+1]][0]:
                    new_sequence[i+1]=new_sequence[i]
                else:
                    new_sequence[i]=new_sequence[i+1]
            elif phenotypes[sequence[i]][0]>phenotypes[sequence[i+1]][0]:
                new_sequence[i+1]=None
                    
        return list(filter(None.__ne__,dict.fromkeys(new_sequence)))
          

    def stripRedundancy(sequences):
        return tuple(set((tuple(seq) for seq in sequences if seq)))

     ##main loop
    for i, (sim,sets) in enumerate(zip(full_simulations,full_sets)):
        sets=[timeCondition(sim,single_set) for single_set in sets]
          
        clean_set=stripRedundancy(sets)
        full_sets[i]=[fset for fset in clean_set if len(fset)>=1]
        if not full_sets[i]:
            full_simulations[i]=None
            full_sets[i]=None
            

def edgeTopology(ep):
    if ep[0] == ep[1]:
        return 1
    elif ep[0]%4 == ep[1]%4:
        return 2
    else:
        return 0

def edgeClassification(ep,h_0):
    topology=edgeTopology(ep)
    if topology==1:
        return 'homodimeric'
    ##true heterodimer
    elif topology==0:
        return 'heterodimeric'
    ##unclear
    else:
        return 'heterodimeric' if h_0 else 'Du-Sp' 

def generateRecord(full_simulations,full_sets,use_raw,samples):
    new_data = defaultdict(list)
    formed_interactions, evolved_interactions = [], []

    
    fail_rate = [0,0,0]
    
    if samples:
        simulation_data = list(zip(full_simulations,full_sets))
        random.shuffle(simulation_data)
        #simulation_data = random.shuffle(list(zip(list(full_simulations),list(full_sets))))
    else:
        simulation_data = zip(full_simulations,full_sets)
        
    for sim, sets in simulation_data:#full_simulations,full_sets):
        trimmed_nodes=[]
        node_details=defaultdict(dict)
        homomeric_discovery, heteromeric_discovery = {}, {}
        heteromeric_compositions = []

        #take longest evolutionary branch to analyse
        branch = sorted(sorted(sets,reverse=True),key=len,reverse=True)[0]
        
        dont_break = False
        fail_rate[0] += 1

        newly_formed_interactions, existing_formed_interactions = [], []
        previous_generation = 0
        initial_details={}

        observed_pairs = set()
        branch = branch[:3]
        for stage, (leaf,parent) in enumerate(zip(branch,(-1,)+branch)):
            
            if len(observed_pairs - {tuple(sorted(e%4 for e in edge[:2])) for edge in sim[leaf][4:]}) >= 1:

                fail_rate[1] += 1
                #print('Unstable evo')
                break

            ##calculate supporting information on homology, occurence, etc.
            minimals=defaultdict(list)
            for edge in sim[leaf][4:]:
                edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                observed_pairs.add(edge_pair)
                minimals[edge_pair].append(edge[2])
                         
            for ep,homol in minimals.items():
                if parent == -1 or ep not in node_details[parent]:
                    node_details[leaf][ep]=[min(homol),sim[leaf][2],0]
                else:
                    node_details[leaf][ep]=node_details[parent][ep][:]
                    node_details[leaf][ep][2]+=1
                                   
            if leaf in trimmed_nodes:
                continue
            else:
                trimmed_nodes.append(leaf)

            ##use supporting information to create data row for this transition
            heteromeric_composition = {}
            if len({tuple(sorted(e%4 for e in edge[:2])) for edge in sim[leaf][4:]}) != stage+1:
                fail_rate[2] += 1
                #print('Doubled topology, bad evo')
                break
            
            for edge in sim[leaf][4:]:
                edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                
                if edge_pair not in heteromeric_composition or heteromeric_composition[edge_pair][1] == 'homodimeric':
                    heteromeric_composition[edge_pair] = (stage,edgeClassification(edge[:2],node_details[leaf][edge_pair][0]))
                              
                ##discovered in wrong order
                if (sim[leaf][2]-node_details[leaf][edge_pair][1])<0:
                    print("negative time of discovery")
                    break
                              
                         
                if edge[0] == edge[1]:
                    if edge_pair not in homomeric_discovery:
                        homomeric_discovery[edge_pair] = (stage,sim[leaf][2]  - use_raw*previous_generation)
                elif edge[0]%4 == edge[1]%4:
                    if edge_pair not in homomeric_discovery:
                        heteromeric_discovery[edge_pair] = (stage,sim[leaf][2]  - use_raw*previous_generation)
                else:
                    if edge_pair not in heteromeric_discovery:
                        heteromeric_discovery[edge_pair] = (stage,sim[leaf][2]  - use_raw*previous_generation)
 
                    
                              

                if node_details[leaf][edge_pair][2] == 0:
                    newly_formed_interactions.append({'stage':stage,'class':edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})
                             
            previous_generation = sim[leaf][2]
            for (stage,composition) in heteromeric_composition.values():
                    existing_formed_interactions.append({'stage':stage,'class':composition})     

        ##after all leafs recorded, take information if valid simulation
        else:
            for (stage,generation) in heteromeric_discovery.values():
                new_data[stage].append(generation)
            for (stage,generation) in homomeric_discovery.values():
                new_data[10+stage].append(generation)

            formed_interactions.extend(newly_formed_interactions)
            evolved_interactions.extend(existing_formed_interactions)

        if samples and fail_rate[0]-(sum(fail_rate[1:]))>=samples:
            break

                         
    print(f'Failed on: unstable ({fail_rate[1]}), double ({fail_rate[2]}), out of {fail_rate[0]} ({(fail_rate[1]+fail_rate[2])/fail_rate[0]:.2f}) [{fail_rate[0]-fail_rate[1]-fail_rate[2]}]')

    return pd.DataFrame(formed_interactions),pd.DataFrame(evolved_interactions),new_data


def makeRecord(S_hat,mu,rate,use_raw,samples):
    sims,sets=readEvolutionRecord(mu,S_hat,rate,FNAME='scripts/')

    cleanRecord(sims,sets)
    
    return generateRecord(filter(None,sims),filter(None,sets),use_raw,samples)

def loadManyRecords(strengths,rates,mu,L,use_raw):
    results = []
     
    for (S_hat, rate) in product(strengths, rates):
        er = EvolutionResult(L,S_hat,mu,rate)
        er.addData(*makeRecord(S_hat, mu, rate,use_raw))
        results.append(er)
          
    if len(results) == 1:
        results = results[0]
    return results

def loadSimulationResults(use_raw=False,samples=None):
    sim_params = (([.83],[0,.05],.00417,60),([.75],[0,.05],.003125,80), ([.74],[0,0.05],.0025,100), ([.7],[0,.05],.0021,120),([.714],[0,.05],.001786,140))
    
    results = []
    for param in sim_params:
        for (S_hat, rate) in product(*param[:2]):
            ER = EvolutionResult(param[3],S_hat,param[2],rate)
            ER.addData(*makeRecord(ER.S_c, ER.mu, ER.dup_rate,use_raw,samples))
            results.append(ER)
        
        #results.extend(loadManyRecords(*d,use_raw=use_raw))
    return results

     
class EvolutionResult(object):
    __slots__ = ('L','S_c','mu','dup_rate','asym','discov_types','composition_types','discov_times','div')
    def __init__(self,L=None,S_c=None,mu=None,dup_rate=None):
        self.L = L
        self.S_c = S_c
        self.mu = mu
        self.dup_rate = dup_rate
        self.asym = formTime(L,S_c)/formTime(L//2,S_c)
         
    def addData(self,discov_types, composition_types, discov_times):
        self.discov_types = discov_types
        self.composition_types = composition_types
        self.discov_times = discov_times
     
    def __repr__(self):
        return f'Evolution result for \'L:{self.L}, S_c:{self.S_c}, Dup: {self.dup_rate}\''



import matplotlib.pyplot as plt

def deWeight(times,init=1):
     heights = list(np.linspace(1/len(times),init,len(times)))
     unique_times = []
     for index,time in enumerate(reversed(times)):
          if time not in unique_times:
               unique_times.append(time)
          else:
               heights.pop(len(times)-index-1)
     return np.array(unique_times[::-1]),np.array(heights)
          
     
from scipy.signal import savgol_filter
def plotDiv2(res,reslice=False,times=(10,)):
     f,(ax,ax2) = plt.subplots(2,sharex=True)
     gg=getGammas()
     for r,r2,c in zip(res[::2],res[1::2],('g','m','r','b','c')):

          
          H_times = np.array(sum((r.discov_times[slice_] for slice_ in times),[]))
          J_times = np.array(sum((r2.discov_times[slice_] for slice_ in times),[]))


          #if r.L == 140 or r.L == 100 or r.L==120:
               #pass
          #     J_times = np.array([x*1.3 for x in r2.discov_times[10]]+r2.discov_times[1]+r2.discov_times[2])
          #if r.L == 60:
          #     H_times = np.array(r.discov_times[10]+[x*1.5 for x in r.discov_times[1]]+[x*1.5 for x in r.discov_times[2]])

          

          H_u,H_h= deWeight(sorted(H_times),len(H_times)/1500)
          J_u,J_h = deWeight(sorted(J_times),len(J_times)/1500)

          LIMIT = 2
          #ax2.plot(H_u[H_h<LIMIT][:-1], savgol_filter(np.ediff1d(H_h[H_h<LIMIT])/np.ediff1d(H_u[H_h<LIMIT]),101,3),c=c)
          #ax2.plot(J_u[J_h<LIMIT][:-1], savgol_filter(np.ediff1d(J_h[J_h<LIMIT])/np.ediff1d(J_u[J_h<LIMIT]),101,3),c=c,ls=':')
          #ax2.plot(H_u[:-1],np.ediff1d(H_h)/np.ediff1d(H_u),c=c,ls='-')
          #ax2.plot(J_u[:-1],np.ediff1d(J_h)/np.ediff1d(J_u),c=c,ls='--')

          
          

          H_spline = InterpolatedUnivariateSpline(H_u,H_h)
          J_spline = InterpolatedUnivariateSpline(J_u,J_h)


          h_space= np.logspace(0,np.log10(H_u[-1]),1000)
          j_space= np.logspace(0,np.log10(J_u[-1]),1000)

          h_spaceF= np.logspace(0,np.log10(H_u[np.argmax(H_h>=2)]),1000)
          j_spaceF= np.logspace(0,np.log10(J_u[np.argmax(J_h>=2)]),1000)
          


          #ax.plot(h_space,H_spline(h_space),c=c,ls='',marker='o',ms=12,markevery=.05)
          #ax.plot(j_space,J_spline(j_space),ls='',c=c,marker='h',ms=12,markevery=.05)
          
          new_start=formTime(r.L//2,r.S_c)/75
          if reslice:
              j_spaceF = j_space[j_space>=new_start]-new_start
          else:
               #ax2.axvline(new_start,c=c)
               j_spaceF = j_space
               
          ax2.plot(j_spaceF,J_spline(j_spaceF)-H_spline(j_spaceF),c=c)


          argmax = np.argmax(J_spline(j_spaceF)-H_spline(j_spaceF))
          ax.plot([j_spaceF[argmax]]*2,[H_spline(j_spaceF[argmax]),J_spline(j_spaceF[argmax])],mfc='none',mec=c,marker='o',ms=12,c=c)
          

          ax.scatter([formTime(r.L//2,r.S_c)/8.25,formTime(r.L,r.S_c)/100/.75],[1,2],c=c)
          #ax.scatter([formTime(r.L//2,r.S_c)/75,formTime(r.L,r.S_c)/200],[1,2],c=c)
          
          ax.plot(H_u,H_h,ls='-',c=c,label=f'{r.L},{r.dup_rate}',alpha=1)#,marker='<',ms=10,markevery=.05)
          
          ax.plot(J_u,J_h,ls=(0,(6,4)),c=c,label=f'{r2.L},{r2.dup_rate}',alpha=1)#,marker='^',ms=10,markevery=.05)

          #ax.plot(np.linspace(0,1e6,101),3e-06*np.linspace(0,1e6,101),'k--')


     #ax.axhline(1,c='k',lw=3,ls='--')
     #ax.axhline(2,c='k',lw=3,ls=':')
          

     #ax.set_yscale('log')
     ax.set_xscale('log')
     ax.legend()

     #ax2.set_yscale('log')
     ax2.set_xscale('log')
     
     plt.show(block=False)

from scipy.interpolate import UnivariateSpline,InterpolatedUnivariateSpline

def getGammas():
     return {60:1.71, 80:3.32, 100: 3.25, 120:4.46, 140:3.74}
