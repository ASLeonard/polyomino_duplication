import sys
#add local paths to load custom methods
if not any(('scripts' in pth for pth in sys.path)):
     sys.path.append('scripts/')
     
from interface_methods import *
import numpy as np

from copy import deepcopy
from sys import argv
from pickle import load,dump
from multiprocessing import Pool
from functools import partial
from collections import defaultdict, Counter
from itertools import combinations, product, groupby
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


from matplotlib.colors import ListedColormap, Normalize, LinearSegmentedColormap
import matplotlib as mpl


def readBinaryVectors(fname,run,L):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/{}_Run{}.txt'.format(fname,run),dtype=np.uint16).reshape(-1,L+1)


def norm_rows(a):
     x=np.sum(a)
     return a if (x==0) else a/x

import matplotlib.colors as mpc
import numpy.ma as ma

def plotNewHomology(run,L,edges):
     
     raw_data = np.fromfile('scripts/{}_Run{}.txt'.format('Zomology',run),dtype=np.uint16).reshape((-1,L+1))
     data = {}
     for i in range(4):
          for j in range(4):
               data[(i,j)] = raw_data[i*4+j::16]

     f,axes= plt.subplots(len(edges),1,sharex=True,sharey=True)

     colors= ['forestgreen','royalblue','orangered']
     for i,edge in enumerate(edges):
          plotHomology(data[edge],ax=axes[i],user_color=colors[i])
          axes[i].axis('off')
          #return data[edge]
     plt.show(block=False)
     
     return
     
 
     
def avg(myArray, N=2):
    cum = np.cumsum(myArray,0)
    result = cum[N-1::N]/float(N)
    result[1:] = result[1:] - result[:-1]

    remainder = myArray.shape[0] % N
    if remainder != 0:
        if remainder < myArray.shape[0]:
            lastAvg = (cum[-1]-cum[-1-remainder])/float(remainder)
        else:
            lastAvg = cum[-1]/float(remainder)
        result = np.vstack([result, lastAvg])

    return result     

def plotHomology(data,L=None, ax=None,user_color='red'):
     if isinstance(data,int):
          data = readBinaryVectors('Zomology',data,L)

     show = ax is None
     if ax is None:
          f,ax = plt.subplots()
          
     data=np.apply_along_axis(norm_rows,1,data.astype(np.float))
     data = avg(data,20)
     pop_grid= ma.masked_equal(data[:300].T,0)

     cm = LinearSegmentedColormap.from_list("", ['gainsboro',user_color])
     px=ax.pcolormesh(pop_grid,cmap=cm,norm=mpc.LogNorm(.005,.25),rasterized=False)
     #plt.colorbar(px,ax=ax)
     #vmin=pop_grid[1:].min(), vmax=pop_grid[1:].max())
     #print(pop_grid[1:].min(),pop_grid[1:].max())
     #
     if show:
          plt.show(block=False)
     


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
          #sets=[growthCondition(sim,single_set) for single_set in sets]
          sets=[timeCondition(sim,single_set) for single_set in sets]
          
        
          clean_set=stripRedundancy(sets)
          full_sets[i]=[fset for fset in clean_set if len(fset)>=1]
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

     new_data = defaultdict(list)
     new_df = []
     new_comp = []
     cnt=0
     terminal_states=[0,0]
     fail_rate = [0,0,0]
     
     for sim, sets in zip(full_simulations,full_sets):
          trimmed_nodes=[]
          node_details=defaultdict(dict)
          heteromeric_discovery = {}
          heteromeric_compositions = []
          sets=[sorted(sorted(sets,reverse=True),key=len,reverse=True)[0]]

          dont_break = False
          fail_rate[0] += 1
          previous_generation = 0
          for branch in sets:
               initial_details={}
               cnt+=1
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
                              continue
                         
                         if stage==0 and True and edge[0] == edge[1]:
                              if edge_pair not in heteromeric_discovery:
                                   heteromeric_discovery[edge_pair] = (len(heteromeric_discovery),sim[leaf][2]  - previous_generation)


                         if node_details[leaf][edge_pair][2] == 0:
                              new_df.append({'stage':stage,'class':edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})
                             
                         DATA.append({'occurence':node_details[leaf][edge_pair][2],'generation':sim[leaf][2],'t_0':sim[leaf][2]-node_details[leaf][edge_pair][1],'homology':edge[2],'edge_pair': edgeTopology(edge[:2]),'h_0':node_details[leaf][edge_pair][0],'class': edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})

  
                    for (stage,composition) in heteromeric_composition.values():
                         new_comp.append({'stage':stage,'class':composition})     
               else:
                    dont_break = True
               previous_generation = sim[leaf][2]
               if not dont_break:
                    #fail_rate[1] += 1
                    break
          ##after all leafs recorded, take information on max phenotype
          else:
               for (stage,generation) in heteromeric_discovery.values():
                    new_data[stage].append(generation)

               if not dont_break:
                    break

               node_details=dict(node_details)

               edges = sim[leaf]
               for edge in edges[4:]:
                    edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                    ##if is a terminal heterodimeric edge, was it initially a homodimeric one
                    if edgeTopology(edge[:2])!=1:
                         terminal_states[node_details[leaf][edge_pair][0] == 0]+=1
                         
     print(f'Failed on: unstable ({fail_rate[1]}), double ({fail_rate[2]}) out of {fail_rate[0]} ({(fail_rate[1]+fail_rate[2])/fail_rate[0]})')
     return pd.DataFrame(DATA),terminal_states,pd.DataFrame(new_df),new_data, pd.DataFrame(new_comp)

def makeRecord(S_hat,mu,rate,dup=True):
     sims,sets=readEvoRecord3(mu,S_hat,rate,dup)

     diversity = calculateDiversity(sims,sets)
     cleanRecord(sims,sets)

     return generateRecord(filter(None,sims),filter(None,sets)), diversity


def loadManyRecords(strengths,rates,mu,L):
     evo_records = {}
     evo_ratios = {}
     evo_diversities = {}
     results = []
     
     for (S_hat, rate) in product(strengths, rates):
          ((evo_records[(S_hat, rate)], *evo_ratios[(S_hat, rate)]), evo_diversities[(S_hat, rate)]) = makeRecord(S_hat, mu, rate)
          er = EvolutionResult(L,S_hat,mu,rate)
          er.addData(evo_ratios[(S_hat, rate)][1],evo_ratios[(S_hat, rate)][3],evo_ratios[(S_hat, rate)][2])
          results.append(er)
          
          
     if len(results) == 1:
          results = results[0]
     return evo_records, evo_ratios, evo_diversities,results

def mergeDataFrames(data_frames):
     new_frames=[df.assign(S_c=S_c,du_rate=du_rate) for (S_c, du_rate), df in data_frames.items()]

     df = pd.concat(new_frames)
     for (S_c, du_rate) in data_frames.items():
          if 'Du-Sp' not in df.loc[(df['S_c']==S_c) & (df['du_rate']==du_rate)]['class']:
               df=df.append({'class':'Du-Sp','edge_pair':0,'generation':0,'h_0':-100,'homology':-100,'occurence':0,'t_0':0,'S_c':S_c,'du_rate':du_rate},ignore_index=True)
               

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
     for  (S_c,du_rate), ratio in bulk_data.items():
          #du_rate=1
          tuples.append((du_rate, S_c, ratio[1]/ratio[0]))
          print(du_rate, S_c, '({})'.format(ratio[1]/ratio[0]),ratio[1],ratio[0])
     vals=np.array(tuples).T
     #print(vals)
     sc=plt.scatter(*vals[:2], c=vals[2],norm=LogNorm(),s=300)
     plt.colorbar(sc)
     #plt.xscale('log')
     plt.show(block=False)
     #return ax,vals

def makeGrid(ax,vs):
     xs=v[0].reshape(3,-1)
     ys=v[1].reshape(3,-1)
     cols=v[2].reshape(3,-1)
     
     CS=plt.contour(xs,ys,cols)
     ax.clabel(CS, fontsize=9, inline=1)
     
def plotTrends(dframes):
     if isinstance(dframes,dict):
          dframes=list(dframes.values())
     elif not isinstance(dframes,list):
          dframes=[dframes]
     f,axes=plt.subplots(len(dframes))
     if len(dframes)==1:
          axes=[axes]
          
     for ax,dframe in zip(axes,dframes):
          sns.pointplot(x="occurence", y="homology", hue="class",capsize=.2, palette="YlGnBu_d", data=dframe,ax=ax,esimator=np.mean,ci=95)
          addRandomExpectations(ax,128)
          
     plt.show(block=False)

def plotTrendsOther(df):

     g = sns.catplot(x='occurence',y='homology',hue='class',data=v2,row='S_c',col='du_rate',kind='point')
     g.set(ylim=(0, None))
     for ax in g.axes.ravel():
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


def doNow(rec):
     rec = np.array(rec[0][2]).T
     return rec
     
     
def makeGriddedAvg(data):
     counters = defaultdict(list)

     for pair in data.T:
          counters[pair[1]].append(pair[0])
     

     results = []
     for k,v in sorted(counters.items()):
          print(k,np.mean(v),len(v))
          results.append([np.mean(v)/(formTime(128,.6875)/250),k])
     plt.scatter(*np.array(results).T)
     
def reconstruct(data):
     times = []
     prev=1
     for row in data:
          if row[1] == 1:
               times.append({(row[1],0):row[0]})
          elif row[1] == prev:
               continue
          else:
               times[-1][(row[1],0)]=row[0]
          prev = row[1]
     return times
               
     

     
def plotDiversity(data_dict):
     f,axes = plt.subplots(2)
     slopes=[]
     for (S_hat, rate), diversity in data_dict.items():
         slopes.append(plotSingleDiversity(axes, diversity, S_hat,rate))

     df = pd.DataFrame(slopes)
     sns.scatterplot(data=df, x='S',y='slope',hue='rate',ax=axes[0])
     axes[0].axhline(0,ls='--',c='darkgray')
     plt.legend()
     axes[0].set_ylabel('Num Phenotype')
     axes[1].set_ylabel('Max size')
     plt.show(block=False)
     
def plotSingleDiversity(axes,data,S,rate):

     rate_lines = {0:':',.0015:'-.',0.015:'-'}
     colours = {K/128:V for (K,V) in zip(range(86,98,2),cm.get_cmap('viridis')(np.linspace(0,1,6)))}
     max_g=max(max(run.values()) for run in data)
     discovs=np.zeros((len(data),max_g,2),dtype=np.uint8)
     for r,simulation in enumerate(data):
          
          hits=sorted([(v,k[0]) for k,v in simulation.items()])
          count_val = 1
          size_val = hits[0][1]
          for start,end in zip(hits,hits[1:]+[(max_g,-1)]):
               discovs[r,start[0]:end[0]]=(count_val,size_val)
               count_val += 1
               size_val =  max(start[1],end[1])

     #f,axes = plt.subplots(2)
     slope_points = []
     def plotHatched(index):
          rescaled_X = np.arange(max_g)/(formTime([64,128][index],S)/100)#getExpectedInteractionTime(32,S)
          
          y=np.mean(discovs[...,index],axis=0)
          y_err=np.std(discovs[...,index],axis=0,ddof=1)
          col=axes[index].plot(rescaled_X,y,lw=2,label=S,ls=rate_lines[rate],color=colours[S])[0].get_color()
          y_low=np.maximum(y-y_err,0)
          slope_points.append({'S':S,'rate':rate,'slope': stats.linregress(rescaled_X[3*y.size//4:], y[3*y.size//4:])[0]})
          #print('{}, {} -> slope: {:.3f}'.format(S,rate, stats.linregress(rescaled_X, y)[0]))
        
          #axes[index].plot(y_low,c=col,ls='--')
          #axes[index].plot(y+y_err,c=col,ls='--')
          #axes[index].fill_between(np.arange(max_g), y_low, y+y_err,alpha=.5,hatch='////',facecolor = 'none',edgecolor=col)

     #rand_samples=np.random.choice(discovs.shape[0], 10, replace=False)
     for index,ls in zip(range(2),('-',':')):

          plotHatched(index)
          #for sample in rand_samples:
          #     axes[index].plot(discovs[sample,...,index],ls=ls,lw=.5)

     
     return slope_points[0]
     

     
          



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


import sys


def plotRidge(df,df2):
     df['L']=0
     df2['L']=1
     df = pd.concat([df,df2])
     sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

     plot_code='occurence'
     df = df.loc[df['class']=='Du-Sp']
     
     # Initialize the FacetGrid object
     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
     g = sns.FacetGrid(df, row=plot_code, hue='L', aspect=10, height=1,sharey=False)#, palette=pal
     
     # Draw the densities in a few steps
     #g.map(sns.kdeplot, "homology", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.1)
     g.map(sns.kdeplot, "homology", clip_on=False, lw=2, bw=.1)
     g.map(plt.axhline, y=0, lw=2, clip_on=False)
     g.map(plt.axvline, x=64, lw=2,color='dimgrey',ls='--',clip_on=False)
     
     
     # Define and use a simple function to label the plot in axes coordinates
     def label(x, color, label):
          ax = plt.gca()
          ax.text(0, .2, label, fontweight="bold", color=color,
                  ha="left", va="center", transform=ax.transAxes)
          

     g.map(label, "occurence")
     # Set the subplots to overlap
     g.fig.subplots_adjust(hspace=-.1)
     
     # Remove axes details that don't play well with overlap
     g.set_titles("")
     g.set(yticks=[])
     g.despine(bottom=True, left=True)
     plt.show(block=False)


def plotImbalance(data1,data2):
     df = []
     tts=['homo','hetero','homo2','hetero2','x','y']
     for j, data in enumerate((data1,data2)):
          for i in range(max(data['stage'])):
               data_S = data.loc[data['stage']==i]
               imbalance = 1#sum(data_S['class']=='homodimeric')+sum(data_S['class']=='heterodimeric')
               df.append({'t':tts[j*2],'stage':i,'imb': (sum(data_S['class']=='homodimeric')/imbalance) if sum(data_S['class']=='homodimeric')> 0 else 0})
               df.append({'t':tts[j*2+1],'stage':i,'imb':-(sum(data_S['class']=='heterodimeric')/imbalance) if sum(data_S['class']=='heterodimeric')> 0 else 0})
               #df.append({'t':tts[j*3+1],'stage':i,'imb': sum(data_S['class']=='Du-Sp')})
                          #'imb':np.log10(imbalance) if imbalance>0 else -np.log10(-imbalance)})
          
     df = pd.DataFrame(df)
     plt.figure()
     sns.barplot(orient='h',data=df,y='stage',hue='t',x='imb')
     plt.show(0)

def plotImbalance2(data1,data2):
     df = []
     tts=['x','homo','hetero','homo2','hetero2']
     for j, data in enumerate((data1,data2)):
          data = data.loc[(data['class']=='Du-Sp') | (data['class']=='heterodimeric')]
          if len(data) == 0:
               continue
          for i in range(max(data['stage'])):
               data_S = data.loc[data['stage']==i]
               if len(data_S) == 0:
                    continue
               imbalance = 1#sum(data_S['class']=='Du-Sp')+sum(data_S['class']=='heterodimeric')

               df.append({'t':tts[j*2],'stage':i,'imb': sum(data_S['class']=='Du-Sp')/imbalance})
               df.append({'t':tts[j*2+1],'stage':i,'imb':-sum(data_S['class']=='heterodimeric')/imbalance})

          
     df = pd.DataFrame(df)
     plt.figure()
     sns.barplot(orient='h',data=df,y='stage',hue='t',x='imb')
     plt.show(0)

def loadRaw(Mu,S,D):
     nums= []
     i=0
     for line in open('/rscratch/asl47/Duplication/EvoRecords/EvoRecord_Mu{:.6f}_S{:.6f}_D{:.6f}.txt'.format(Mu,S,D)):
          if i==0:
               nums.append(int(line.split()[2]))

          i+=1

          i = i%3
     return nums

from interface_formation import formTime
def plotCumin(data,ls='-'):
     total_runs = len(data[0])
     #plt.figure()
     for stage in data.keys():
          if stage>2:
               continue
          times = np.array([1]+sorted(data[stage]))
          plt.plot(times/(formTime(128,.71875)/100),np.arange(0,len(times))/total_runs,label=stage,ls=ls)

     plt.xscale('log')
     plt.legend()
     plt.show(0)

#from itertools import cycle
from scipy.stats import expon,gamma

def plotTimex(*datas):
     N = 1
     f,ax = plt.subplots(2)
     pop = 100

     cmap = cm.get_cmap('copper')

     

     asyms = []
     combined_data = []
     for data in datas:
     #     for i in data.discov_times.keys():
     #               continue
     #          data.discov_times[0].extend(data.discov_times[i])
          asyms.append(formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c))

     scaler=mpc.LogNorm(min(asyms),max(asyms))

     for stage in range(N):
          L_scaler = 2-stage
          for data in datas:
               #sns.kdeplot(np.array(data.discov_times[stage])/formTime(data.L//L_scaler,data.S_c)*pop,cut=0,ax=ax[stage],ls='-' if data.dup_rate else '-.',label=f'L:{data.L}, S_c:{data.S_c}, Dup:{data.dup_rate}',color=cmap(scaler(formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c))),kernel='epa')
               print(len(np.array(data.discov_times[stage])/formTime(data.L//L_scaler,data.S_c)*pop))
               fit_p = expon.fit(np.array(data.discov_times[stage])/formTime(data.L//L_scaler,data.S_c)*pop)
               print(data.L,fit_p)
               ax[stage].plot(np.linspace(0,10,101),expon(*fit_p).pdf(np.linspace(0,10,101)),marker='h',markevery=.1,label=f'L:{data.L}, S_c:{data.S_c}, Dup:{data.dup_rate}',c=cmap(scaler(formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c))))
               #
               sns.distplot(np.array(data.discov_times[stage])/formTime(data.L//L_scaler,data.S_c)*pop,bins=np.linspace(.001,10,10),ax=ax[stage],color=cmap(scaler(formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c))),kde=False,hist_kws={'histtype':'step','density':1})
          


          ax[stage].set_yscale('log',nonposy='clip')
          #ax[stage].set_xscale('log',nonposx='clip')
     ax[0].legend()
     plt.show(0)

#from matplotlib.sankey import Sankey
def chartSankey(data):
    for stage in range(max(data.discov_types['stage'])+1):
        print('Incoming at stage ',stage)
        for t_class in ('homodimeric','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.discov_types['class']==t_class) & (data.discov_types['stage']==stage)))

          
        print('Compositions')
        for t_class in ('homodimeric','Du-Sp','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.composition_types['class']==t_class) & (data.composition_types['stage']==stage)))
          

class EvolutionResult(object):
    __slots__ = ('L','S_c','mu','dup_rate','asym','discov_types','composition_types','discov_times')
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

def getRes():
     loadD = (([.83],[0,.05],.00417,60),([.75],[0,.05],.003125,80), ([.74],[0,0.05],.0025,100), ([.7],[0,.05],.0021,120),([.714],[0,.05],.001786,140))
     res=[]
     for d in loadD:
          er, err, ed, (QQ)= loadManyRecords(*d)
          res.extend(QQ)
     return res
