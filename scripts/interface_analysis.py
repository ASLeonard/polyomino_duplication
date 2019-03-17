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

def plotBs(a):
     plt.figure()
     c=['r','b','g']
     for i in a:
          for g,j in enumerate(i.T):
               plt.plot(range(1000),j,c[g])
     plt.show(block=False)

from matplotlib.colors import ListedColormap, Normalize
import matplotlib as mpl

def homm(run,L):
     return [[int(i) for i in line.split()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Yomology_Run{}.txt'.format(run))]
def strr(run,L):
     return [[L-int(i) for i in line.split()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Strengths_Run{}.txt'.format(run))]

def readBinaryVectors(fname,run,L):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/{}_Run{}.txt'.format(fname,run),dtype=np.uint16).reshape(-1,L+1)


def norm_rows(a):
     x=np.sum(a)
     return a if (x==0) else a/x

import matplotlib.colors as mpc
import numpy.ma as ma
def plotHom(run,L,norm=True,annotate=False):

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
          fixed_pids=[]#{tuple(i) for i in getFixedPhenotypes(LoadPIDHistory(run))}
          for pid,details in annotations.items():
               alph=1 if pid in fixed_pids else 1
                    
               for edge in details[2:]:
                    ax.scatter(details[0],edge[2],c=[cm.tab20((edge[0]%4*4+edge[1]%4)/16)],alpha=alph)
     plt.show(block=False)

def readEvoRecord(run):
     lines=[line.rstrip() for line in open('/scratch/asl47/Data_Runs/Bulk_Data/Evo_Run{}.txt'.format(run))]
     d={}
     for l in lines:
          parts=l.split()
          d[tuple(int(i) for i in parts[:2])]=tuple(int(i) for i in parts[2:4])+tuple([tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)])
     return d
def readEvoRecord3(mu,S_c,rate,duplicate=True):
     lines=[line.rstrip() for line in open('/rscratch/asl47/Duplication/EvoRecords/EvoRecord_Mu{:.6f}_S{:.6f}_{}{:.6f}.txt'.format(mu,S_c,'D' if duplicate else 'I',rate))]
     print(lines)
     simulations=[[]]
     sets=None
     for line in lines:
          if line == '':
               print(sets)
               print("new run")
               
          elif ',' in line:
               sets=getSuperSets([{int(i) for i in bm.split()} for bm in line.split(',') if bm])
               
          else:
               parts=line.split()
               #key=tuple(int(i) for i in parts[:2])
               simulations[-1].append(tuple(int(i) for i in parts[:4])+tuple([tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)]))
               
     return simulations,sets

def generateRecord(sim,sets):
     orderedOcc={}
     
def readEvoRecord2(mu,S_c,rate,duplicate=True):
     lines=[line.rstrip() for line in open('/rscratch/asl47/Duplication/EvoRecords/EvoRecord_Mu{:.6f}_S{:.6f}_{}{:.6f}.txt'.format(mu,S_c,'D' if duplicate else 'I',rate))]
     d=[{}]
     #print(lines)
     prev_phen_size=(0,0)
     for l in lines:
          parts=l.split()
          key=tuple(int(i) for i in parts[:2])
          if key <= prev_phen_size:
               d.append({})
          d[-1][key]=tuple(int(i) for i in parts[2:4])+tuple([tuple(int(i) for i in parts[q:q+4])+(float(parts[q+4]),) for q in range(4,len(parts)-4,5)])
          
          prev_phen_size=key
               
               
     return [dict_run for dict_run in d if dict_run]

from scipy import stats
import pandas as pd
def compileEvoRecord(runs):
     
     data=[]
     cnt=-1
     for simulation in runs:
          cnt+=1
          occurences=defaultdict(int)
          initial_homologies={}
          initial_edges={}
          dp={v[0]:v[2:] for v in sorted(simulation.values())}
          for gen,v in dp.items():
               minimals=defaultdict(list)
               for edge in v:
                    edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                    minimals[edge_pair].append(edge[2])
                    
               for ep,homol in minimals.items():
                    if ep not in initial_homologies:
                         initial_homologies[ep]=min(homol)
                    if ep not in initial_edges:
                         initial_edges[ep]=gen
               
               for edge in v:
                    edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                    if occurences[edge_pair] == 0 and edge[2] != initial_homologies[edge_pair]:
                         if cnt in []:
                              continue
                         #break
                    #if occurences[edge_pair]>6:
                    #     print(simulation)
                         #print(cnt)
                         #break
                         #print(edge_pair,occurences[edge_pair], edge[2], initial_homologies[edge_pair])
                         #return simulation
                         

                    data.append({'occurence':occurences[edge_pair],'generation':gen,'t_0':gen-initial_edges[edge_pair],'homology':edge[2],'edge_pair':edge_pair[0]==edge_pair[1],'h_0':initial_homologies[edge_pair]})


               for co in minimals.keys():
                    occurences[co]+=1
     return pd.DataFrame(data)

def readCompiledEvoRecord(add_type='Dup', T='High'):
     return pd.read_csv('/rscratch/asl47/EvoRecords/EvoRecord_{}_{}.pd'.format(add_type,T))

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
          f.colorbar(px,ax=ax)#ax.colorbar(px)
     plt.show(block=False)
     
def plotEvoRecord(df,key='occurence'):
     f,ax=plt.subplots()

     for homologues in (0,1):
          d=df.loc[df['edge_pair']==homologues]
          plt.scatter(d[key],d['homology'],c=d['occurence'],alpha=0.5,cmap='plasma',marker=('o','s')[homologues])#,norm=Normalize(0,90))

     homol_data=defaultdict(list)
     rand_data=defaultdict(list)
     for xv,h_0,hom in zip(df[key],df['h_0'],df['homology']):
          if h_0<=40:
               #if xv!=0 and hom==0:
               #     continue
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
     #slope, intercept, r_value, p_value, std_err=stats.linregress(sorted(homol_data),[np.mean(homol_data[k]) for k in sorted(homol_data)])
     #plt.plot(range(min(homol_data),max(homol_data)+1),[slope*x+intercept for x in range(min(homol_data),max(homol_data)+1)],'r--')

     print("ratio is {}/{}".format(len(df.loc[(df['occurence']==0) & (df['h_0']<=40)]),len(df.loc[(df['occurence']==0) & (df['h_0']>40)])))
     ax.set_ylabel('homology',fontsize=24)
     ax.set_xlabel(key,fontsize=24)
     plt.show(block=False)

def getSuperSets(list_of_sets):
     super_sets=[]
     for l1 in range(len(list_of_sets)):
          for l2 in range(l1+1,len(list_of_sets)):
               if list_of_sets[l1] < list_of_sets[l2]:
                    break
          else:
               super_sets.append(list_of_sets[l1])
     return super_sets
     
def plotDiversity(data,max_g):
     discovs=np.zeros((len(data),max_g,2),dtype=np.uint8)
     for r,simulation in enumerate(data):
          
          hits=sorted([(v[0],k[0]) for k,v in simulation.items()])
          count_val = 1
          try:
               size_val = hits[0][1]
          except:
               print(hits,r)
          for start,end in zip(hits,hits[1:]+[(max_g,-1)]):
               discovs[r,start[0]:end[0]]=(count_val,size_val)
               count_val += 1
               size_val =  max(start[1],end[1])
     #return discovs
     #plt.figure()
     def plotHatched(index):
          y=np.mean(discovs[...,index],axis=0)
          y_err=np.std(discovs[...,index],axis=0,ddof=1)
          col=plt.plot(y,lw=2,label=id(data))[0].get_color()
          y_low=np.maximum(y-y_err,0)
        
          #plt.plot(y_low,c=col,ls='--')
          #plt.plot(y+y_err,c=col,ls='--')
          #plt.fill_between(range(max_g), y_low, y+y_err,alpha=.5,hatch='////',facecolor = 'none',edgecolor=col)

     
     for index,ls in zip(range(2),('-',':')):
          plotHatched(index)
          #plt.fill_between(range(max_g), y-y_err, y+y_err,alpha=.2)
          #plt.plot(,ls=ls)
     plt.legend()
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
     #mutations=np.loadtxt('/scratch/asl47/Data_Runs/Bulk_Data/Mutation_Run{}.txt'.format(run),dtype=np.uint8)[:selections.shape[0]-1,:].reshape(-1)
     #homologies=np.genfromtxt('/scratch/asl47/Data_Runs/Bulk_Data/Homology_Run1.txt',dtype=np.float64,delimiter=",")

     sizes=np.loadtxt('/scratch/asl47/Data_Runs/Bulk_Data/Size_Run{}.txt'.format(run),dtype=np.uint8).reshape(-1)

     

     cols=np.array(['k','darkgreen','darkred','blue','gainsboro'])
     lws=np.array([0.5,1,1,1])
     lc = mc.LineCollection(lines, linewidths=1,linestyle='-',cmap='inferno')
     
     if colorize=='M':
          lc.set_array(mutations)
     elif colorize=='H':
          lc.set_array(homologies)
     elif colorize=='S':
          lc.set_array(sizes)     
     
     plt.colorbar(lc)
     ax.add_collection(lc)

         
from itertools import chain     
def plotPhen(pids_raw,run,thresh=0.25):
     pids=pids_raw

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

if __name__ == "__main__":
     main()

def main():
     args=sys.argv
     compileEvoRecord()
