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

from scipy import stats
import pandas as pd

#GLOBAL PIDS
null_pid,init_pid=np.array([0,0],dtype=np.uint8),np.array([1,0],dtype=np.uint8)


from matplotlib.colors import ListedColormap, Normalize
import matplotlib as mpl


def readBinaryVectors(fname,run,L):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/{}_Run{}.txt'.format(fname,run),dtype=np.uint16).reshape(-1,L+1)


def norm_rows(a):
     x=np.sum(a)
     return a if (x==0) else a/x

import matplotlib.colors as mpc
import numpy.ma as ma
def plotHomologyEvolution(run,L,norm=True,annotate=False):

     f,axes=plt.subplots(2,1,sharex=True)
     for ax,func in zip(axes,('Zomology','Strengths')):

          data=readBinaryVectors(func,run,L)
          if norm:
               data=np.apply_along_axis(norm_rows,1,data.astype(np.float))
          pop_grid= ma.masked_equal(data.T,0)

          px=ax.pcolormesh(pop_grid,cmap='RdGy',norm=mpc.LogNorm(vmin=pop_grid.min(), vmax=pop_grid.max()))

     axes[0].set_ylabel('Homology')
     axes[1].set_ylabel('Strength')
     
     f.colorbar(px,ax=axes)
     if annotate:
          annotations=readEvoRecord(run)
          ax=axes[0]
          fixed_pids={tuple(i) for i in getFixedPhenotypes(LoadPIDHistory(run))}
          for pid,details in annotations.items():
               alph=1 if pid in fixed_pids else 1
                    
               for edge in details[2:]:
                    ax.scatter(details[0],edge[2],c=[cm.tab20((edge[0]%4*4+edge[1]%4)/16)],alpha=alph)
     plt.show(block=False)


def readEvoRecord3(mu,S_c,rate,duplicate=True):
     lines=[line.rstrip() for line in open('/rscratch/asl47/Duplication/EvoRecords/EvoRecord_Mu{:.6f}_S{:.6f}_{}{:.6f}.txt'.format(mu,S_c,'D' if duplicate else 'I',rate))]


     simulations=[[]]
     sets=[]
     for line in lines:
          if line == '':
               if simulations[-1]:
                    simulations.append([])
               
          elif ',' in line:
               sets.append(getSuperSets([{int(i) for i in bm.split()} for bm in line.split(',') if bm]))
               
          else:
               parts=line.split()
               simulations[-1].append(tuple(int(i) for i in parts[:4])+tuple(tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)))
               
     return simulations[:-1],sets

def cleanRecord(full_simulations,full_sets):

     ##REDUNDANT 
     def growthCondition(phenotypes, sequence):
          if not sequence:
               return None
          if len(sequence) == 1:
               return sequence
          elif phenotypes[sequence[-1]][0]<=phenotypes[sequence[-2]][0]:
               del sequence[-1]
               growthCondition(phenotypes,sequence)
          else:
               return sequence

     def timeCondition(phenotypes, sequence):
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
          
          #{item for sublist in sequences for item in sublist}
          return tuple(set((tuple(seq) for seq in sequences if seq)))
     #def stripTrivials(sequences):
     #     return [i for i in a if len(i)>1]tuple(set((tuple(seq) for seq in sequences)))

     ##main loop
     for i, (sim,sets) in enumerate(zip(full_simulations,full_sets)):
          #sets=[growthCondition(sim,single_set) for single_set in sets]
          sets=[timeCondition(sim,single_set) for single_set in sets]
          
        
          clean_set=stripRedundancy(sets)
          full_sets[i]=[fset for fset in clean_set if len(fset)>1]
          if not full_sets[i]:
               full_simulations[i]=None
               full_sets[i]=None
                        
def generateRecord(full_simulations,full_sets):

     def edgeTopology(ep):
          if ep[0] == ep[1]:
               return 1
          elif ep[0]%4 == ep[1]%4:
               return 2
          else:
               return 0

     def edgeClassification(ep,h_0):
          ##true homodimer
          topology=edgeTopology(ep)
          if topology==1:
               return 'homodimeric'
          ##true heterodimer
          elif topology==0:
               return 'heterodimeric'
          ##unclear
          else:
               return 'heterodimeric' if h_0 else 'Du-Sp' 
          
     DATA=[]
     #cnt=0
     terminal_states=[0,0]
     for sim, sets in zip(full_simulations,full_sets):
          trimmed_nodes=[]
          node_details=defaultdict(dict)
          for branch in sets:
               initial_details={}
               for leaf,parent in zip(branch,(-1,)+branch):

                    ##calculate supporting information on homology, occurence, etc.
                    minimals=defaultdict(list)
                    for edge in sim[leaf][4:]:
                         edge_pair=tuple(sorted(e%4 for e in edge[:2]))
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
                    for edge in sim[leaf][4:]:
                         edge_pair=tuple(sorted(e%4 for e in edge[:2]))

                         ##discovered in wrong order
                         if (sim[leaf][2]-node_details[leaf][edge_pair][1])<0:
                              print("negative time of discovery")
                              continue

                         DATA.append({'occurence':node_details[leaf][edge_pair][2],'generation':sim[leaf][2],'t_0':sim[leaf][2]-node_details[leaf][edge_pair][1],'homology':edge[2],'edge_pair': edgeTopology(edge[:2]),'h_0':node_details[leaf][edge_pair][0],'class': edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})
                         
          ##after all leafs recorded, take information on max phenotype
          else:
               node_details=dict(node_details)
               #print(cnt)
               #cnt+=1
               leaf_index, edges = max(enumerate(sim), key=itemgetter(1))
               for edge in edges[4:]:
                    edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                    ##if is a terminal heterodimeric edge, was it initially a homodimeric one
                    if edgeTopology(edge[:2])!=1:
                         terminal_states[node_details[leaf_index][edge_pair][0] == 0]+=1
                        
               
     return pd.DataFrame(DATA),terminal_states

def makeRecord(S_hat,mu,rate,dup=True):
     sims,sets=readEvoRecord3(mu,S_hat,rate,dup)
     cleanRecord(sims,sets)

     return generateRecord(filter(None,sims),filter(None,sets)), calculateDiversity(filter(None,sims),filter(None,sets))

def calculateDiversity(sims,sets):
     diversity=[]
     for sim, set_ in zip(sims,sets):
          phen_discovery={}
          for pathway in set_:
               for step in pathway:
                    pid=sim[step][:2]
                    if pid in phen_discovery:
                         phen_discovery[pid] = min(phen_discovery[pid], sim[step][2])
                    elif not phen_discovery or pid[0] >= max(phen_discovery)[0]:
                         phen_discovery[pid] = sim[step][2]
                         
          diversity.append(phen_discovery)
     return diversity

def plotE2(df,norm=True):
     f,axes=plt.subplots(2,1)
     data=np.zeros((2,2500,129))
     
     for row in df.itertuples(index=False):
          data[int(row[0]),row[1],row[3]]+=1

     if norm:
          data=np.apply_along_axis(norm_rows,1,data.astype(np.float))
     pop_grid= ma.masked_equal(data,0)

     for index,ax in enumerate(axes):
          px=ax.pcolormesh(pop_grid[index].T,cmap='RdGy',norm=mpc.LogNorm(vmin=pop_grid[index].min(), vmax=pop_grid[index].max()))
          f.colorbar(px,ax=ax)
     plt.show(block=False)
     
def plotEvoRecord(df,key='occurence'):
     f,ax=plt.subplots()

     for homologues in (0,1,2):
          d=df.loc[df['edge_pair']==homologues]
          plt.scatter(d[key],d['homology'],c=d['h_0'],alpha=0.5,cmap='cividis',marker=('d','P','.')[homologues],norm=Normalize(0,100))

     homol_data=defaultdict(list)
     rand_data=defaultdict(list)
     for xv,h_0,hom in zip(df[key],df['h_0'],df['homology']):
          if h_0<=40:
               if xv!=0 and hom==0:
                    continue
               homol_data[xv].append(hom)
          else:
               rand_data[xv].append(hom)

     def running_mean(x, N):
          cumsum = np.cumsum(np.insert(x, 0, 0)) 
          return (cumsum[N:] - cumsum[:-N]) / float(N)

     window=20 if key!='occurence' else 2
     for group,c in zip((rand_data,homol_data),('k','r')): 
          plt.plot(sorted(group)[max(0,window//2-1):-window//2],running_mean([np.median(group[k]) for k in sorted(group)],window),c=c,lw=2)

     if key == 't_0':
          plt.plot(64*(1-np.exp(-np.arange(max(df['t_0']))/2/(256))),ls='--',c='c')

     print("ratio is {}/{}".format(len(df.loc[(df['occurence']==0) & (df['h_0']<=40)]),len(df.loc[(df['occurence']==0) & (df['h_0']>40)])))
     ax.set_ylabel('homology',fontsize=24)
     ax.set_xlabel(key,fontsize=24)
     plt.show(block=False)

     #sns.stripplot(x='occurence',y='homology',hue='class',data=dframe)
     #sns.violinplot(x='occurence',y='homology',hue='class',data=dframe)

def plotPhaseSpace(bulk_data):
     f,ax=plt.subplots()
     tuples=[]
     for (du_rate, S_c), (_,ratio) in bulk_data.items():
          tuples.append((du_rate, S_c, ratio[1]/ratio[0]))
          print(du_rate, S_c, '({})'.format(ratio[1]/ratio[0]),ratio[1],ratio[0])
     vals=np.array(tuples).T
     #print(vals)
     sc=plt.scatter(*vals[:2], c=vals[2],norm=LogNorm(),s=300)
     plt.colorbar(sc)
     plt.xscale('log')
     plt.show(block=False)
     return ax,vals

def makeGrid(ax,vs):
     xs=v[0].reshape(3,-1)
     ys=v[1].reshape(3,-1)
     cols=v[2].reshape(3,-1)
     
     CS=plt.contour(xs,ys,cols)
     ax.clabel(CS, fontsize=9, inline=1)
     
def plotTrends(dframes):
     if not isinstance(dframes,list):
          dframes=[dframes]
     f,axes=plt.subplots(len(dframes))
     if len(dframes)==1:
          axes=[axes]
          
     for ax,dframe in zip(axes,dframes):
          sns.pointplot(x="occurence", y="homology", hue="class",capsize=.2, palette="YlGnBu_d", data=dframe,ax=ax,esimator=np.mean,ci=95)
          addRandomExpectations(ax,128)
          
     plt.show(block=False)
     
def addRandomExpectations(ax,L_I):
     ax.axhline(L_I/2,c='#666666',ls='--',lw=2,zorder=0)
     for sign in (-1,1):
          ax.axhline(L_I/2+sign*np.sqrt(L_I)/2,c='#666666',ls='--',lw=.75,zorder=0)

def getSuperSets(list_of_sets):
     super_sets=[]
     for l1 in range(len(list_of_sets)):
          for l2 in range(l1+1,len(list_of_sets)):
               if list_of_sets[l1] < list_of_sets[l2]:
                    break
          else:
               super_sets.append(sorted(list_of_sets[l1]))
     return super_sets
     
def plotDiversity(data):
     max_g=max(max(run.values()) for run in data)
     discovs=np.zeros((len(data),max_g,2),dtype=np.uint8)
     for r,simulation in enumerate(data):
          
          hits=sorted([(v,k[0]) for k,v in simulation.items()])
          count_val = 1
          try:
               size_val = hits[0][1]
          except:
               print(hits,r)
          for start,end in zip(hits,hits[1:]+[(max_g,-1)]):
               discovs[r,start[0]:end[0]]=(count_val,size_val)
               count_val += 1
               size_val =  max(start[1],end[1])

     f,axes = plt.subplots(2)
     def plotHatched(index):
          y=np.mean(discovs[...,index],axis=0)
          y_err=np.std(discovs[...,index],axis=0,ddof=1)
          col=axes[index].plot(y,lw=2,label=id(data))[0].get_color()
          y_low=np.maximum(y-y_err,0)
        
          axes[index].plot(y_low,c=col,ls='--')
          axes[index].plot(y+y_err,c=col,ls='--')
          axes[index].fill_between(range(max_g), y_low, y+y_err,alpha=.5,hatch='////',facecolor = 'none',edgecolor=col)

     rand_samples=np.random.choice(discovs.shape[0], 10, replace=False)
     for index,ls in zip(range(2),('-',':')):
          plotHatched(index)
          for sample in rand_samples:
               axes[index].plot(discovs[sample,...,index],ls=ls,lw=.5)

     #plt.legend()
     axes[0].set_ylabel('Num Phenotype')
     axes[1].set_ylabel('Max size')
     plt.show(block=False)

     
          



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

     #plt.figure()
     for pid,details in data.items():
          for edge in details[2:]:
               plt.scatter(edge[2],edge[3],c=[cm.tab20((edge[0]%4*4+edge[1]%4)/16)])
     
     plt.show(block=False)
                                        
     
def lazy(run,style='S'):
     add_selection_layer(plotPhen2(LoadPIDHistory(run)),LSHB(run,100),run,style)
     
def Int(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Interactions_Run{}.txt'.format(run))][g].split('.')[c]
def Bint(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Binteractions_Run{}.txt'.format(run))][g].split(',')[c]
def Hom(run,g,c):
     return [line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Homology_Run{}.txt'.format(run))][g].split(',')[c]





         

     
from scipy.stats import binom
def getExpectedInteractionTime(L_I,S_c):
     binomial=binom(L_I,.5)
     thresh=math.floor(S_c*L_I)-1
     return 1/binomial.sf(thresh)
##baseline is .001 mutation

import subprocess
import sys
def runnerCore(S_CRITICAL,DU_FLAG='J'):
     runs={.671875:250,.6875:350,.703125:500, .71875:625, .734375:750,.75:1000}

     default_args=' -N 10 -P 100 -B 100 -X .25 -F 1 -A 1 -V 0 -T 10 -Y {} -M .001 -D {} -G {}'.format(S_CRITICAL, runs[S_CRITICAL], runs[S_CRITICAL]*10)

     dead_flag={'J':'K','K':'J'}

     for G_RATE in (.0001,.001,.01):
          special_args=' -L {} -{} {} -{} 0'.format(G_RATE, DU_FLAG, G_RATE, dead_flag[DU_FLAG])
          subprocess.run('../bin/DuplicationEvolution -E'+default_args+special_args,shell=True)     
     
     
def main():
     if len(sys.argv)==2:
          runnerCore(float(sys.argv[1]))
     print('Done okay')
     
if __name__ == "__main__":
     main()
