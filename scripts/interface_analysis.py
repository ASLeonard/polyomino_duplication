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
import math
import glob

import warnings

import matplotlib.pyplot as plt
import matplotlib.cm as cm

#GLOBAL PIDS
null_pid,init_pid=np.array([0,0],dtype=np.uint8),np.array([1,0],dtype=np.uint8)

def parallelAnalysis(S_star,t,mu,gamma,runs,offset=0,run_code='F'):
     setBasePath('scratch')
     chosen_function=analysePhylogenetics if run_code=='F' else analyseHomogeneousPopulation
     pool = Pool()
     data_struct=pool.map(partial(chosen_function, S_star,t,mu,gamma), range(offset,offset+runs)) 
     pool.close()
     if run_code=='F':
          dump(data_struct, open('Y{}T{}Mu{}J{}K{}L{}O{}.pkl'.format(mu,S_star,t,gamma,offset), 'wb'))
     else:
          np.savez_compressed('Mu{}Y{}T{}F{}O{}'.format(mu,S_star,t,gamma,offset),data_struct)
          


def collateAnalysis(params,runs):
     N_samps=10
     full_data=[]
     for r in runs:
          try:
               full_data.append(load(open('Y{}T{}Mu{}J{}K{}L{}O{}.pkl'.format(*params+(r,)), 'rb')))
          except:
               print('missing pickle for run ',r)
     N_runs=0
     full_transition=defaultdict(lambda: defaultdict(int))

     for phen_tran in filter(None,full_data):
          N_runs+=1
          for (phen_in,phen_out),count in phen_tran.items():
               full_transition[phen_in][phen_out]+=count
          
     return (N_runs,convertDoubleNestedDict(full_transition))
     
def convertRaggedArray(list_of_lists):
     long_length=sorted([len(lis) for lis in list_of_lists],reverse=True)[len(list_of_lists)//2]
     rect_arr=np.empty((len(list_of_lists),long_length))
     for i,strs in enumerate(list_of_lists):
          rect_arr[i]=strs[:long_length]+[np.nan]*(long_length-len(strs[:long_length]))
     return rect_arr

def analysePhylogenetics(run,params,full_pIDs=False):
     #s,p,st,phen_table=LoadAll(run,params)
     s=LSHB(run,params)
     p=LPB(run,params)
     phen_table=LoadPhenotypeTable(run)
     phen_table[(2,2)]=(2,1,1,9)
     transitions=KAG(p,s)
     if not transitions:
          print("Empty run at {}".format(run))
          return None
          
     if full_pIDs:
          transitions={(phen_table[k[0]],phen_table[k[1]]):cnt for k,cnt in transitions.items()}
     return transitions
     transitions=ret_val[1]
     failed_transitions=ret_val[2]
     bond_data=treeBondStrengths(ret_val[0],st)
     if full_pIDs:
          transitions={(phen_table[k[0]],phen_table[k[1]]):cnt for k,cnt in transitions.items()}
          failed_transitions={(phen_table[k[0]],phen_table[k[1]]):cnt for k,cnt in failed_transitions.items()}
          bond_data={phen_table[k]:v for k,v in bond_data.items()}
     return (bond_data,transitions,failed_transitions)

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
          
def KAG(phenotypes_in,selections):
     phenotypes=phenotypes_in.copy()
     max_gen,pop_size=selections.shape
     
     forest,temp_forest=[],[]
     transitions=defaultdict(int)

     def __growDescendentTree(tree,max_depth=float('inf')):
          gen_val=tree.gen
          descendents=tree.seq[0]
          while gen_val<(max_gen-1):
               
               new_descendents=[]
               for descendent in descendents:
                    new_descendents.extend([child for child in np.where(selections[gen_val]==descendent)[0] if np.array_equal(phenotypes_in[gen_val+1,child],tree.pID)])

               if math.isinf(max_depth):
                    phenotypes[gen_val,descendents]=null_pid
                                   
               if not new_descendents:
                    break
               if (gen_val-tree.gen)>=max_depth:
                    return True
               descendents=new_descendents
               tree.seq.append(descendents)
               gen_val+=1
          else:
               phenotypes[gen_val,descendents]=null_pid
                              
     def __addBranch():

          pid_ref=phenotypes_in[g_idx,c_idx]

          transitions[tuple(tuple(_) for _ in (pid_ref,init_pid if g_idx==0 else phenotypes_in[g_idx-1,p_idx]))]+=1

          temp_forest.append((True,Tree(pid_ref,0,0,g_idx,[[c_idx]])))

          return True 
               
     for C_INDEX in range(pop_size):
          if np.array_equal(phenotypes[max_gen-1,C_INDEX],null_pid):
               continue
          
          c_idx=C_INDEX
          g_idx=max_gen-1
          p_idx=selections[g_idx-1,c_idx]
          pid_ref=phenotypes[max_gen-1,c_idx]
         
          
          while g_idx>0:
               if np.array_equal(phenotypes[g_idx-1,p_idx],null_pid):
                    if np.array_equal(phenotypes_in[g_idx-1,p_idx],pid_ref):
                         temp_forest.append((False,Tree(pid_ref,0,(-1,-1),g_idx,[[c_idx]])))
                    else:
                         if not __addBranch():
                              return None
                         break
               
               elif not np.array_equal(phenotypes[g_idx-1,p_idx],pid_ref):
                    if not __addBranch():
                         return None
                    pid_ref=phenotypes[g_idx-1,p_idx]
                    
          
               g_idx-=1
               c_idx=p_idx
               p_idx=selections[g_idx-1,p_idx]
          else:
               if not np.array_equal(pid_ref,init_pid) and not __addBranch():
                    return None
          
          while temp_forest:
               (alive,tree)=temp_forest.pop()
               __growDescendentTree(tree)
               if alive:
                    forest.append(tree)  

     return dict(transitions)
     
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
          bond_data[tuple(tree.pID)].append((tree.gen,dict(bond_maps)))
     return dict(bond_data)     

def plotBs(a):
     plt.figure()
     c=['r','b','g']
     for i in a:
          for g,j in enumerate(i.T):
               plt.plot(range(1000),j,c[g])
     plt.show(block=False)

from matplotlib.colors import ListedColormap
import matplotlib as mpl

def homm(run,L):
     return [[int(i) for i in line.split()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Yomology_Run{}.txt'.format(run))]

def homm2(fname,run,L):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/{}_Run{}.txt'.format(fname,run),dtype=np.uint16).reshape(-1,L+1)

def strr(run,L):
     return [[L-int(i) for i in line.split()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Strengths_Run{}.txt'.format(run))]


def norm_rows(a):
     x=np.sum(a)
     return a if (x==0) else a/x

import matplotlib.colors as mpc
import numpy.ma as ma
def plotHom(run,L,norm=True,anno=False):

     f,axes=plt.subplots(2,1,sharex=True)
     for ax,func in zip(axes,('Zomology','Strengths')):
          if func==strr:
               continue
               print(homm)
          data=homm2(func,run,L)
          if norm:
               data=np.apply_along_axis(norm_rows,1,data.astype(np.float))
          pop_grid= ma.masked_equal(data.T,0)
          
          #pop_grid=np.ma.zeros((L+1,len(data)))
          #pop_grid.mask=True
          #for i,row in enumerate(data):
          #     c=Counter(row)c
          #     for k,v in c.items():
          #          pop_grid.mask[k,i]=False
          #          pop_grid[k,i]=v/sum(c.values())
          px=ax.pcolormesh(pop_grid,cmap='RdGy',norm=mpc.LogNorm(vmin=pop_grid.min(), vmax=pop_grid.max()))

     axes[0].set_ylabel('Homology')
     axes[1].set_ylabel('Strength')
     
     f.colorbar(px,ax=axes)
     if anno:
          annotes=lls(run)
          ax=axes[0]
          fixed_pids={tuple(i) for i in getFixedPhenotypes(LoadPIDHistory(run))}
          for pid,details in annotes.items():
               alph=1 if pid in fixed_pids else .2
                    
               for edge in details[2:]:
                    ax.scatter(details[0],edge[2],c=[cm.tab20((edge[0]%4*4+edge[1]%4)/16)],alpha=alph)
     plt.show(block=False)

def lls(run):
     lines=[line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Evo_Run{}.txt'.format(run))]
     d={}
     for l in lines:
          parts=l.split()
          
          d[tuple(int(i) for i in parts[:2])]=tuple(int(i) for i in parts[2:4])+tuple([tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)])
     return d

def getFixedPhenotypes(pids):
     def perGeneration(pid_slice):
          cc=Counter(list(chain.from_iterable(pid_slice.flat)))
          for bad_pid in ((0,0),(255,0)):
               if bad_pid in cc:
                    del cc[bad_pid]
          return max(cc, key=lambda key: cc[key])
     return np.apply_along_axis(perGeneration,1,pids)

def scatL(run):
     data=lls(run)

     plt.figure()
     for pid,details in data.items():
          for edge in details[2:]:
               plt.scatter(edge[2],edge[3],c=[cm.tab20((edge[0]%4*4+edge[1]%4)/16)])
     
     plt.show(block=False)
                                        
     
def lazy(run,style='S'):
     add_selection_layer(plotPhen2(LoadPIDHistory(run)),LSHB(run,250),run,style)
     
def Int(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Interactions_Run{}.txt'.format(run))][g].split('.')[c]
def Bint(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Binteractions_Run{}.txt'.format(run))][g].split(',')[c]
def Hom(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Homology_Run{}.txt'.format(run))][g].split(',')[c]

def plotPhen2(pids_raw):
     pids=ObjArray([[tuple(i) for i in row] for row in pids_raw])
     ref_pids=list(np.unique(pids))

     dets=[pid for pid in ref_pids if (0,0) not in pid]
     non_dets=[pid for pid in ref_pids if (0,0) in pid]

     if ((255,0),) in dets:
          non_dets.append(dets.pop(dets.index(((255,0),))))
         
     cmap_base = cm.get_cmap('tab20b',len(dets))
     cmap_base2 = cm.get_cmap('tab20c',20)
     cmap_rare=[]

     counter=0
     for rare in non_dets:
          raw_rare=tuple(sub_pid for sub_pid in rare if sub_pid!=(0,0))
          if raw_rare in dets:
               index = dets.index(raw_rare)
               cmap_rare.append(cmap_base.colors[index].copy())
               cmap_rare[-1][3]=.25
          elif len(raw_rare)==0:
               cmap_rare.append((1,1,1,1))
          elif raw_rare==((255,0),):
               cmap_rare.append((0,0,0,1))
          else:
               cmap_rare.append(cmap_base2.colors[counter%20].copy())
               cmap_rare[-1][3]=.25
               counter+=1


     cmap_full = ListedColormap(np.vstack((cmap_base.colors,cmap_rare)))
     #print(cmap_full.colors)
     c={K:i for i,K in enumerate(dets+non_dets)}
     #print(c)
     
     pop_grid=np.empty(pids.shape).T
     for i,j in np.ndindex(pids.shape):
          pop_grid[j,i]=c[pids[i,j]]
          
     #len(ref_pids))

     
     f,ax=plt.subplots()
     #img = plt.imshow(np.array([list(range(len(dets)))]), cmap=cmap_base)
     #img.set_visible(False)

     #cbar=plt.colorbar(drawedges=True)#orientation="vertical")

     px=plt.pcolormesh(pop_grid,cmap=cmap_full)
     
     sm=cm.ScalarMappable(cmap=cmap_base)
     sm.set_array(np.linspace(0,1,len(dets)))
     
     cbar=plt.colorbar(sm,drawedges=True)
     #tick_locs = (np.arange(len(ref_pids)) + 0.5)*(len(ref_pids)-1)/len(ref_pids)
     tick_locs = (np.arange(len(dets)) + 0.5)/len(dets)
     cbar.set_ticks(tick_locs)

     #cbar.ax.set_yticklabels(dets+non_dets, fontsize=14, weight='bold')
     cbar.ax.set_yticklabels(dets, fontsize=14, weight='bold')
     plt.show(block=False)
     return ax

from matplotlib import collections  as mc
def add_selection_layer(ax,selections,run,colorize='M'):
     lines=[[(-.5,i+.5,),(.5,i+.5)] for i in range(selections.shape[1])]
     for g_ind, p_ind in np.ndindex(selections.shape):
          if g_ind==(selections.shape[0]-1):
               continue
          lines.append([(g_ind+.5,selections[g_ind,p_ind]+.5),(g_ind+1.5,p_ind+.5)])
     mutations=np.loadtxt('/scratch/asl47/Data_Runs/Bulk_Data/Mutation_Run{}.txt'.format(run),dtype=np.uint8)[:selections.shape[0]-1,:].reshape(-1)
     #homologies=np.genfromtxt('/scratch/asl47/Data_Runs/Bulk_Data/Homology_Run1.txt',dtype=np.float64,delimiter=",")

     sizes=np.loadtxt('/scratch/asl47/Data_Runs/Bulk_Data/Size_Run{}.txt'.format(run),dtype=np.uint8).reshape(-1)

     ##[line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Interactions_Run1.txt')][815].split('.')[96]##
     

     cols=np.array(['k','darkgreen','darkred','blue','gainsboro'])
     lws=np.array([0.5,1,1,1])
     lc = mc.LineCollection(lines, linewidths=1,linestyle='-',cmap='inferno')#color=col_opts)#lws[mutations]
     
     if colorize=='M':
          lc.set_array(mutations)
     elif colorize=='H':
          lc.set_array(homologies)
     elif colorize=='S':
          lc.set_array(sizes)     
     
     plt.colorbar(lc)
     ax.add_collection(lc)

def add_duplication_layer(ax):
     dups=[[int(i) for i in line.rstrip()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Homology_Run0.txt')]
     

     
from itertools import chain     
def plotPhen(pids_raw,run,thresh=0.25):
     pids=pids_raw#ObjArray([[tuple(i) for i in row] for row in pids_raw])

     unique_pids=[(0,0),(1,0),(255,0)]+list(LoadPhenotypeTable(run).keys())
     gens,pop_size=pids_raw.shape

     gen_counts={K:[0]*gens for K in unique_pids}
     for i,row in enumerate(pids):
          for pid in chain.from_iterable(row.flat):
               gen_counts[pid][i]+=1

     
     plt.figure()
     
     for k,v in gen_counts.items():
          if max(v)<pop_size*thresh:
               continue
          plt.plot(range(gens),v,label=k,ls='-' if k!=(0,0) else ':',c=None if k!=(0,0) else 'k')
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

     
               
def main(argv):
     model_type=int(argv[2])
     if argv[1]=='internal':
          
          HPC_FLAG=argv[3]=='1'
          run=int(argv[4])
          format_params=tuple(float(i) for i in argv[5:11])
          run_params=int(argv[11])# if HPC_FLAG else format_params
          
          if model_type==1 or model_type==0:
               with open('Y{}T{}Mu{}J{}K{}L{}O{}.pkl'.format(*format_params+(run,)),'wb') as f:
                    dump(analysePhylogenetics(run,run_params,1),f)
               for used_file in glob.glob('*Run{}*'.format(run)):
                    pass
                    #os.remove(used_file)
                    
          else:
               print("hi")

     elif argv[1]=='external':
          format_params=tuple(float(i) for i in argv[3:9])
          file_pth='/rscratch/asl47/Pickles/Y{}T{}Mu{}J{}K{}L{}'.format(*format_params)
          run_gen=range(int(argv[9]))
          if model_type==1 or model_type==0:
               with open(file_pth+'.pkl', 'wb') as f:
                    dump(collateAnalysis(format_params,runs=run_gen), f)
          else:
               print("hi")
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

