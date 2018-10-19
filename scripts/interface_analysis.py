from interface_methods import *
import numpy as np

from copy import deepcopy
from sys import argv
from pickle import load,dump
from multiprocessing import Pool
from functools import partial
from collections import defaultdict,Counter
from itertools import combinations,product,groupby
from operator import itemgetter

def parallelAnalysis(S_star,t,mu,gamma,runs,offset=0,run_code='F'):
     setBasePath('scratch')
     chosen_function=analysePhylogenetics if run_code=='F' else analyseHomogeneousPopulation
     pool = Pool()
     data_struct=pool.map(partial(chosen_function, S_star,t,mu,gamma), range(offset,offset+runs)) 
     pool.close()
     if run_code=='F':
          dump(data_struct, open('Mu{}Y{}T{}F{}O{}.pkl'.format(mu,S_star,t,gamma,offset), 'wb'))
     else:
          np.savez_compressed('Mu{}Y{}T{}F{}O{}'.format(mu,S_star,t,gamma,offset),data_struct)
          
def collateNPZs(S_star,t,mu,gamma,runs):
     full_data=[]
     for r in runs:
          full_data.extend(np.load(open('Mu{}Y{}T{}F{}O{}.npz'.format(mu,S_star,t,gamma,r), 'rb'))['arr_0'])
     return np.asarray(full_data)

def collateAnalysis(S_star,t,mu,gamma,runs):
     N_samps=10
     full_data=[]
     for r in runs:
          try:
               full_data.append(load(open('Mu{}Y{}T{}F{}O{}.pkl'.format(mu,S_star,t,gamma,r), 'rb')))
          except:
               print('missing pickle for run ',r)
     N_runs=0
     full_transition=defaultdict(lambda: defaultdict(int))
     failed_jumps=defaultdict(lambda: defaultdict(int))
     raw_evolutions=defaultdict(lambda: defaultdict(list))
     average_evolutions=defaultdict(dict)
     sample_evolutions=defaultdict(dict)
     
     for str_evo,phen_tran,fails in filter(None,full_data):
          N_runs+=1
          for (phen_in,phen_out),count in phen_tran.items():
               full_transition[phen_in][phen_out]+=count
          for (phen_in,phen_out),count in fails.items():
               failed_jumps[phen_in][phen_out]+=count
          for phen,evo_list in str_evo.items():
               for (gen,evos) in evo_list:
                    for (bond_key,new_bond,_,_),strs in evos.items():
                         raw_evolutions[phen][(bond_key,new_bond)].append(strs)
          
     for phen in raw_evolutions.keys():
          samps=None
          for bond_key,evos in raw_evolutions[phen].items():
               str_arr=convertRaggedArray(evos)
               if samps is None or max(samps)>=str_arr.shape[0]:
                    samps=np.random.choice(str_arr.shape[0],min(str_arr.shape[0],N_samps),replace=False)
               average_evolutions[phen][bond_key]=np.nanmean(str_arr,axis=0)
               sample_evolutions[phen][bond_key]=str_arr[samps]
     return (N_runs,)+tuple(convertDoubleNestedDict(tr) for tr in (full_transition,failed_jumps)) + tuple(dict(res) for res in (average_evolutions,sample_evolutions))
     
def convertRaggedArray(list_of_lists):
     long_length=sorted([len(lis) for lis in list_of_lists],reverse=True)[len(list_of_lists)//2]
     rect_arr=np.empty((len(list_of_lists),long_length))
     for i,strs in enumerate(list_of_lists):
          rect_arr[i]=strs[:long_length]+[np.nan]*(long_length-len(strs[:long_length]))
     return rect_arr

def analysePhylogenetics(S_star,t,mu,gamma,r):
     s,p,st,phen_table=LoadAll(mu=mu,t=t,run=r,S_star=S_star,gamma=gamma)
     forest,transitions,failed_jumps=KAG(p,s,st)
     if not forest:
          return None
     bond_data=treeBondStrengths(forest,st)
     return (bond_data,transitions,failed_jumps)

class Tree(object):
     __slots__ = ('pID','bonds','new_bond','gen','seq')
     def __init__(self,pid=None,bonds=None,new_bond=None,gen=None,seq=None):
          self.pID=pid
          self.bonds=bonds
          self.new_bond=new_bond
          self.gen=gen
          self.seq=seq
     def __repr__(self):
          return '{},{}'.format(self.pID,self.gen)
          
def KAG(phenotypes_in,selections,interactions):
     phenotypes=phenotypes_in.copy()
     max_gen,pop_size=phenotypes.shape
     
     null_pid,init_pid=(0,0),(1,0)
     BAD_return=([],{},{})
     forest,temp_forest=[],[]
     transitions=defaultdict(int)

     def __growDescendentTree(tree,max_depth=float('inf')):
          gen_val=tree.gen
          descendents=tree.seq[0]
          while gen_val<(max_gen-1):
               new_descendents=[]
               for descendent in descendents:
                    new_descendents.extend([child for child in np.where(selections[gen_val]==descendent)[0] if (phenotypes_in[gen_val+1,child]==tree.pID and interactions[gen_val+1,child].bonds==tree.bonds)])
                    phenotypes[gen_val,descendent]=null_pid
                    
               if not new_descendents:
                    break
               if (gen_val-tree.gen)>=max_depth:
                    return True
               descendents=new_descendents
               tree.seq.append(descendents)
               gen_val+=1
          else:
               for d in descendents:
                    phenotypes[gen_val,d]=null_pid
     def __followTree(gen_val,c_idx,duration):
          descendents=[]
          for _ in range(duration):
               new_descendents=[]
               for descendent in descendents:
                    new_descendents.extend([child for child in np.where(selections[gen_val]==descendent)[0] if phenotypes_in[gen_val+1,child]==tree.pID])
               if not new_descendents:
                    break
               descendents=new_descendents
               tree.seq.append(descendents)
               gen_val+=1
               
                    
     def __addBranch():
          bond_ref=interactions[g_idx,c_idx].bonds
          pid_ref=phenotypes_in[g_idx,c_idx]
          if g_idx:
               new_bond=list(set(bond_ref)-set(interactions[g_idx-1,p_idx].bonds))
               transitions[(pid_ref,phenotypes_in[g_idx-1,p_idx])]+=1
          else:
               new_bond=bond_ref
               transitions[(pid_ref,init_pid)]+=1
          temp_forest.append((True,Tree(pid_ref,bond_ref,new_bond[0],g_idx,[[c_idx]])))
          return True if len(new_bond)==1 else False
               
     for C_INDEX in range(pop_size):
          if phenotypes[max_gen-1,C_INDEX]==null_pid:
               continue
          
          c_idx=C_INDEX
          g_idx=max_gen-1
          p_idx=selections[g_idx-1,c_idx]
          pid_ref=phenotypes[max_gen-1,c_idx]
          bond_ref=interactions[g_idx,c_idx].bonds
          
          while g_idx>0:
               if phenotypes[g_idx-1,p_idx]==null_pid:
                    if phenotypes_in[g_idx-1,p_idx]==pid_ref:
                         temp_forest.append((False,Tree(pid_ref,bond_ref,(-1,-1),g_idx,[[c_idx]])))
                    else:
                         if not __addBranch():
                              return BAD_return
                         break
               
               elif phenotypes[g_idx-1,p_idx]!=pid_ref:
                    if not __addBranch():
                         return BAD_return
                    bond_ref=interactions[g_idx-1,p_idx].bonds
                    pid_ref=phenotypes[g_idx-1,p_idx]
               elif interactions[g_idx-1,p_idx].bonds !=bond_ref:
                    temp_forest.append((False,Tree(pid_ref,bond_ref,(-1,-1),g_idx,[[c_idx]])))
                    bond_ref=interactions[g_idx-1,p_idx].bonds
          
               g_idx-=1
               c_idx=p_idx
               p_idx=selections[g_idx-1,p_idx]
          else:
               if pid_ref!=init_pid and not __addBranch():
                    return BAD_return
          
          while temp_forest:
               (alive,tree)=temp_forest.pop()
               __growDescendentTree(tree)
               if alive:
                    forest.append(tree)

     phen_priority={(0,0):0,(1,0):1,(2,0):2,(4,0):2,(4,1):3,(8,0):3,(12,0):4,(16,0):4}
     max_phens=np.vectorize(phen_priority.__getitem__)(phenotypes_in)
     
     #max_priority=np.max(max_phens,axis=1)
     failed_jumps=defaultdict(int)
     for g_idx,c_idx in product(range(max_gen-2,-1,-1),range(pop_size)):
          pid_c=phenotypes[g_idx,c_idx]
          pid_d=(phenotypes_in[g_idx-1,selections[g_idx-1,c_idx]] if g_idx>0 else init_pid)
          if pid_c!=null_pid and pid_c!=pid_d:
              if np.count_nonzero(max_phens[g_idx]>=phen_priority[pid_c])<(pop_size//10) and not __growDescendentTree(Tree(pid_c,interactions[g_idx,c_idx].bonds,(-1,-1),g_idx,[[c_idx]]),2):
                   failed_jumps[(pid_c,pid_d)]+=1
     return (forest,dict(transitions),dict(failed_jumps))
     
def treeBondStrengths(KAG,interactions):
     bond_data=defaultdict(list)
     for tree in KAG: 
          bond_maps=defaultdict(list)
          max_pop=0
          for generation,populations in enumerate(tree.seq,tree.gen):
               if len(populations)<(max_pop//10) and max_pop>(interactions.shape[1]//10):
                    #print(len(populations),max_pop)
                    break
               max_pop=max(max_pop,len(populations))
               inner_bond_maps=defaultdict(list)
               for species in populations:
                    all_bonds=interactions[generation,species].bonds
                    new_bond_type=getBondType(tree.new_bond,all_bonds)
                    for bond,strength in interactions[generation,species]:
                         inner_bond_maps[(getBondType(bond,all_bonds),new_bond_type)+bond].append(strength)

               for k,v in inner_bond_maps.items():
                    bond_maps[k].append(np.mean(v))
          bond_data[tree.pID].append((tree.gen,dict(bond_maps)))
     return dict(bond_data)     

         
def plotPhen2(ps):
     plt.figure()
     refs=list(np.unique(ps))
     c={K:i for i,K in enumerate(refs)}
     z=np.zeros(ps.shape)
     for i,j in product(range(ps.shape[0]),range(ps.shape[1])):
          z[i,j]=c[ps[i,j]]
     plt.pcolormesh(z.T)
     plt.colorbar()
     plt.show(block=False)
     
def plotPhen(ps):
     plt.figure()
     refs=list(np.unique(ps))
     d={K:[] for K in refs}
     for row in ps:
          c=Counter(row)
          for k,v in c.items():
               d[k].append(v)
          for j in refs:
               if j not in c:
                    d[j].append(0)

     for k,v in d.items():
          plt.plot(range(len(list(d.values())[0])),d[k],label=k)
     plt.legend()
     #plt.yscale('log',nonposy='mask')
     plt.show(block=False)


def convertDoubleNestedDict(dict_in):
     return {k:dict(v) for k,v in dict_in.items()}

def getBondType(bond,bonds):
     if checkBranchingPoint(bond,bonds):
          if bond[0]//4==bond[1]//4:
               return 4 #same tile BP
          else:
               return 3 #external BP
     else:
          if bond[0]//4==bond[1]//4:
               return 2 #internal loop
          else:
               return 1 #external SIF

def checkBranchingPoint(bond,bonds):
     test_bonds=list(sum((b for b in bonds if b!= bond), ()))
     return any(x in test_bonds for x in bond)

def consecutiveRanges(data):
     return [map(itemgetter(1), g) for _,g in groupby(enumerate(data), lambda kv:(kv[0]-kv[1]))]

def allUniqueBonds(bonds):
     seen = set()
     return not any(i in seen or seen.add(i) for i in list(sum(bonds, ())))

def writeIt():
     np.savez_compressed('/rscratch/asl47/Pickles/test.npy',[a,b])
def analyseHomogeneousPopulation(S_star,t,mu,gamma,r):
     s,phenotypes,st,phen_table=LoadAll(mu=mu,t=t,run=r,S_star=S_star,gamma=gamma)
     
     max_gen,pop_size=phenotypes.shape
     param_trajectory=[]
     for generation in range(max_gen):
          params=defaultdict(list)
          for species in range(pop_size):
               if phenotypes[generation,species]==(0,0):
                    continue
               if len(st[generation,species].bonds)!=3:
                    continue
               try:
                    bonds={getBondType(bond,st[generation,species].bonds):strength for bond,strength in st[generation,species]}

                    params['a'].append((bonds[4]/bonds[2])**t)
                    params['b'].append((bonds[2]/bonds[3])**t)
               except:
                    return None
          param_trajectory.append([np.mean(params['a']),np.mean(params['b'])])
     return np.asarray(param_trajectory)
     
               
def main(argv):
     if argv[1]=='reduced':
          params=tuple(float(i) for i in argv[2:6])+(int(argv[6]),)
          print("running main for ",argv[6])
          dump(analysePhylogenetics(*params), open('Mu{2}Y{0}T{1}F{3}O{4}.pkl'.format(*params), 'wb'))
          #parallelAnalysis(*(float(i) for i in argv[2:6]),runs=int(argv[6]),offset=int(argv[7]),run_code='F')
     elif argv[1]=='collate':
          with open('/rscratch/asl47/Pickles/Y{}T{}Mu{}F{}.pkl'.format(*(float(i) for i in argv[2:6])), 'wb') as f:
               dump(collateAnalysis(*(float(i) for i in argv[2:6]),runs=range(int(argv[6]))), f)
               
     elif argv[1]=='final':          
          parallelAnalysis(*(float(i) for i in argv[2:6]),runs=int(argv[6]),offset=int(argv[7]),run_code='B')
     elif argv[1]=='group':
          np.savez_compressed('/rscratch/asl47/Pickles/Y{}T{}Mu{}F{}'.format(*(float(i) for i in argv[2:6])),collateNPZs(*(float(i) for i in argv[2:6]),runs=range(0,int(argv[6]),int(argv[7]))))
               
     else:
          print('unknown')
                       
     return

if __name__ == '__main__':
    main(argv)

def PhenotypicTransitions(phen_trans,N=40,crit_factor=0.5):
     common_transitions=deepcopy(phen_trans)
     for phen_key,trans in phen_trans.items():
          for tran,count in trans.items():
               if count<N*crit_factor:
                    del common_transitions[phen_key][tran]

     for key in common_transitions.keys():
          if not common_transitions[key]:
               del common_transitions[key]
     return common_transitions

