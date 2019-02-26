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
          fixed_pids=[]#{tuple(i) for i in getFixedPhenotypes(LoadPIDHistory(run))}
          for pid,details in annotes.items():
               alph=1 if pid in fixed_pids else 1
                    
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

def orderedOcc(runs):
     f,ax=plt.subplots()
     for r in runs:
          occ=defaultdict(list)
          dp={v[0]:v[2:] for v in lls(r).values()}
          for v in dp.values():
               for edge in v:
                    occ[tuple(e%4 for e in edge[:2])].append(edge[2])

          for oc in occ.values():
               c='r' if oc[0]<40 else 'k'
               if len(oc)>10:
                    print(r)
               if len(oc)>1:
                    ax.plot(oc,marker='o',c=c)
               else:
                    ax.scatter([0],oc,c=c)
     plt.show(block=False)
     
def orderedOcc2(runs):
     f,ax=plt.subplots()
     points=[]
     for r in runs:
          for k,v in lls(r).items():
               for edge in v[2:]:
                    points.append((v[0],edge[2]))
     d=np.asarray(points)
     plt.scatter(*d.T)

#dp={v[0]:v[2:] for v in lls(r).values()}
#          for v in dp.values():
 #              for edge in v:
  #                  occ[tuple(e%4 for e in edge[:2])].append(edge[2])

   #       for oc in occ.values():
    #           c='r' if oc[0]<40 else 'k'
     #          if len(oc)>10:
      #              print(r)
       #        if len(oc)>1:
        #            ax.plot(oc,marker='o',c=c)
         #      else:
          #          ax.scatter([0],oc,c=c)
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

