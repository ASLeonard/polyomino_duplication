import sys
#add local paths to load custom methods
if not any(('scripts' in pth for pth in sys.path)):
     sys.path.append('scripts/')
     
from interface_methods import *

import numpy as np
import numpy.ma as ma
from scipy.stats import binom, expon, anderson_ksamp, linregress

from sys import argv
from pickle import load,dump

from collections import defaultdict, Counter
from itertools import product

from scipy.interpolate import UnivariateSpline,InterpolatedUnivariateSpline

import math


import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, Normalize, LogNorm, LinearSegmentedColormap

import pandas as pd

def plotNewHomology(run,L,edges):
    raw_data = np.fromfile(f'Homology_Run{run}.txt',dtype=np.uint16).reshape((-1,L+1))
    data = {}
    for I in range(16):
        data[(I//4,I%4)] = raw_data[I::16]
            
    f,axes= plt.subplots(len(edges),1,sharex=True,sharey=True)
    colors= ['forestgreen','royalblue','orangered','k','k']
    for i,edge in enumerate(edges):
        plotSingleHomology(data[edge],L=L,ax=axes[i],user_color=colors[i])

    plt.show(block=False)
 
def plotSingleHomology(data,L,ax,user_color='red',smooth_window=0):

    data = np.apply_along_axis(lambda a: a if np.sum(a)==0 else a/np.sum(a),1,data.astype(np.float))
    if smooth_window:
        cumulative_smooth = np.cumsum(data,0)
        data = cumulative_smooth[smooth_window-1::smooth_window]/smooth_window
        data = data[1:] - data[:-1]

    ##reverse homology so 100% is at top of plot
    pop_grid= ma.masked_equal(data[:300,::-1].T,0)

    cm = LinearSegmentedColormap.from_list("", ['gainsboro',user_color])
    ax.pcolormesh(pop_grid,cmap=cm,norm=LogNorm(.005,.25),rasterized=False)

def readEvoRecord3(mu,S_c,rate,duplicate=True,FNAME='/rscratch/asl47/Duplication/EvoRecordsz/'):

    def getSuperSets(list_of_sets):
        super_sets=[]
        for l1 in range(len(list_of_sets)):
            for l2 in range(l1+1,len(list_of_sets)):
                if list_of_sets[l1] < list_of_sets[l2]:
                    break
            else:
                super_sets.append(sorted(list_of_sets[l1]))
        return super_sets
    
    lines=[line.rstrip() for line in open('{}EvoRecord_Mu{:.6f}_S{:.6f}_{}{:.6f}.txt'.format(FNAME,mu,S_c,'D' if duplicate else 'I',rate))]

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
use_raw = False          
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
          homomeric_discovery, heteromeric_discovery = {}, {}
          heteromeric_compositions = []
          sets=[sorted(sorted(sets,reverse=True),key=len,reverse=True)[0]]

          dont_break = False
          fail_rate[0] += 1
          
          for branch in sets:
               previous_generation = 0
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
                              break
                              continue
                         
                         if edge[0] != edge[1]:
                              if stage >= 1 and edge_pair not in heteromeric_discovery:
                                   heteromeric_discovery[edge_pair] = (stage,sim[leaf][2]  - use_raw*previous_generation)
 

                         else:
                              if stage == 0 and edge_pair not in homomeric_discovery:
                                   homomeric_discovery[edge_pair] = (stage,sim[leaf][2]  - use_raw*previous_generation)
                              


                         if node_details[leaf][edge_pair][2] == 0:
                              new_df.append({'stage':stage,'class':edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})
                             
                         DATA.append({'occurence':node_details[leaf][edge_pair][2],'generation':sim[leaf][2],'t_0':sim[leaf][2]-node_details[leaf][edge_pair][1],'homology':edge[2],'edge_pair': edgeTopology(edge[:2]),'h_0':node_details[leaf][edge_pair][0],'class': edgeClassification(edge[:2],node_details[leaf][edge_pair][0])})

                    previous_generation = sim[leaf][2]
                    for (stage,composition) in heteromeric_composition.values():
                         new_comp.append({'stage':stage,'class':composition})     
               else:
                    dont_break = True
               
               if not dont_break:
                    #fail_rate[1] += 1
                    break
          ##after all leafs recorded, take information on max phenotype
          else:
               for (stage,generation) in heteromeric_discovery.values():
                    new_data[stage].append(generation)
               for (stage,generation) in homomeric_discovery.values():
                    new_data[10+stage].append(generation)

               if not dont_break:
                    break

               node_details=dict(node_details)

               edges = sim[leaf]
               for edge in edges[4:]:
                    edge_pair=tuple(sorted(e%4 for e in edge[:2]))
                    ##if is a terminal heterodimeric edge, was it initially a homodimeric one
                    if edgeTopology(edge[:2])!=1:
                         terminal_states[node_details[leaf][edge_pair][0] == 0]+=1
                         
     print(f'Failed on: unstable ({fail_rate[1]}), double ({fail_rate[2]}) out of {fail_rate[0]} ({(fail_rate[1]+fail_rate[2])/fail_rate[0]}) [{fail_rate[0]-fail_rate[1]-fail_rate[2]}]')
     return pd.DataFrame(DATA),terminal_states,pd.DataFrame(new_df),new_data, pd.DataFrame(new_comp)

def makeRecord(S_hat,mu,rate,dup=True):
     sims,sets=readEvoRecord3(mu,S_hat,rate,dup,FNAME='scripts/')

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
          er.addData(evo_ratios[(S_hat, rate)][1],evo_ratios[(S_hat, rate)][3],evo_ratios[(S_hat, rate)][2],evo_diversities[(S_hat,rate)])
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
          px=ax.pcolormesh(pop_grid[index].T,cmap='RdGy',norm=LogNorm(vmin=pop_grid[index].min(), vmax=pop_grid[index].max()))
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
     for data in data_dict:#(S_hat, rate), diversity in data_dict.items():
         diversity = data.div
         S_hat= data.S_c
         rate = data.dup_rate
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
          slope_points.append({'S':S,'rate':rate,'slope': linregress(rescaled_X[3*y.size//4:], y[3*y.size//4:])[0]})
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





         

     

def getExpectedInteractionTime(L_I,S_c):
     binomial=binom(L_I,.5)
     thresh=math.floor(S_c*L_I)-1
     return 1/binomial.sf(thresh)
##baseline is .001 mutation






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


def plotTimex(*datas,fit_func=expon,renormalise=True,full_renorm=False,row2=False):
     
     c_thresh=0.75

     np.random.seed(17528175)
     if full_renorm:
          assert renormalise, 'conflicting settings'
     N = 2
     f,ax = plt.subplots(N)
     pop = 100

     

     
     markers={60:'o',80:'s',100:'d',120:'X',140:'*'}
     asyms = []
          
     for data in datas:
          
          asyms.append(formTime(data.L,data.S_c)/formTime(data.L//2,data.S_c))#*getGammas()[data.L]/100)
     print(asyms)

     scaler=LogNorm(min(asyms),max(asyms)*2)

     for stage in range(N):
          mid_data= []
          for asym_val, data in zip(asyms,datas):
               L_scaler = 2-stage
               ##all homomers
               if stage == 0:
                    mutate_rate_adjust = .33 #4 
                    combinatoric_adjust = 4
               ##heteromers
               else:
                    mutate_rate_adjust = 3/4 #3/4 # 5.33
                    combinatoric_adjust = 1#4 #4/36
                         
               gamma_factor = 1
               if full_renorm:
                    if stage == 1 and data.dup_rate > .01:
                         #mutate_rate_adjust = 3/4
                         gamma_factor = getGammas()[data.L]/100
                         #combinatoric_adjust = 1 #4/6
                         L_scaler = 2

               pop = 100 
               raw_data = np.array(data.discov_times[10 if stage == 0 else 1]+([] if (stage == 0 or not row2) else data.discov_times[2]))
               raw_data = np.random.choice(raw_data,size=1500)

               renorm_factor = (formTime(data.L//L_scaler,data.S_c) / pop / mutate_rate_adjust / gamma_factor * combinatoric_adjust) if renormalise else 1

               data_scaled = raw_data/renorm_factor
               if (data.L == 140 or data.L == 100 or data.L==120) and stage == 0 and data.dup_rate>0:
                    data_scaled *=1.3
               if (data.L == 60) and stage == 1 and data.dup_rate==0:
                    data_scaled *=1.5



                    
               cmap = cm.get_cmap('Blues_r' if data.dup_rate>0 else 'Oranges_r')
               asym_color = cmap(scaler(asym_val))
               
               if fit_func:
                    fit_p = fit_func.fit(data_scaled,floc=0)
                    #print(fit_p)
                    x_points = np.linspace(0,max(data_scaled),300)
                    ax[stage].plot(x_points,fit_func(*fit_p).pdf(x_points),marker=markers[data.L],markevery=[0,-1],c=asym_color,ls='-' if data.dup_rate ==0 else '--',lw=4,mew=5,mfc='none',ms=30,alpha=1,label=f'L:{data.L}, S_c:{data.S_c}' if data.dup_rate !=-1 else None)
               if renormalise:
                    if not full_renorm and stage == 1 and data.dup_rate > 0.01:
                         BINS = np.linspace(0,.5,11)
                    else:
                         BINS = np.linspace(0,3,31) if stage ==1 else np.linspace(0,3,31)
               else:
                    BINS=30
                    
               #sns.distplot(data_scaled,bins=BINS,ax=ax[stage],kde=False,hist_kws={'cumulative':True,'histtype':'step','density':1,'lw':2,'alpha':.8,'ls':'-' if data.dup_rate ==0 else '-.'},color=asym_color,label=f'L:{data.L}, S_c:{data.S_c}' if data.dup_rate ==0 else None)#,label=f'{data.dup_rate},{data.S_c}'
               mid_data.append(data_scaled)


          ax[stage].set_yscale('log',nonposy='clip')
          #print(anderson_ksamp(mid_data))
     ax[0].legend()
     ax[0].tick_params(axis='both', which='major', labelsize=16)
     ax[1].tick_params(axis='both', which='major', labelsize=16)
     plt.show(0)


#from matplotlib.sankey import Sankey

def chartSankey(data):
    for stage in range(max(data.discov_types['stage'])+1):
        print('Incoming at stage ',stage)

        for t_class in ('homodimeric','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.discov_types['class']==t_class) & (data.discov_types['stage']==stage)))#/len(data.composition_types.loc[data.composition_types['stage']==stage])*100)
          
        print('Compositions')
        for t_class in ('homodimeric','Du-Sp','heterodimeric'):
            print(f'\t{t_class}: ',sum((data.composition_types['class']==t_class) & (data.composition_types['stage']==stage)))#/len(data.composition_types.loc[data.composition_types['stage']==stage])*100)
    return [sum((data.composition_types['class']==t_class) & (data.composition_types['stage']==2)) for t_class in ('Du-Sp','heterodimeric')]
            
def stackedBars():
     XS=np.arange(3)
     plt.figure()
     new_homomers = np.array((97,3,.5))
     new_heteromers = np.array((3,97,99.5))
     old_homomers = np.array((0,97,100))
     old_heteromers = np.array((0,3,100))
     plt.bar(XS,new_homomers)
     plt.bar(XS,new_heteromers,bottom=new_homomers+old_heteromers+old_homomers)
     plt.bar(XS,old_heteromers,bottom=new_homomers+old_homomers)
     plt.bar(XS,old_homomers,bottom=new_homomers)

     plt.figure()
     new_homomers = np.array((97,87,60))
     new_heteromers = np.array((3,13,40))
     old_homomers = np.array((0,4,12))
     old_dusp = np.array((0,93,93+79))
     old_heteromers = np.array((0,3,16))
     plt.bar(XS,new_homomers)
     plt.bar(XS,new_heteromers,bottom=new_homomers+old_dusp+old_heteromers+old_homomers)
     plt.bar(XS,old_dusp,bottom=new_homomers+old_homomers,hatch='//',color='none')
     plt.bar(XS,old_heteromers,bottom=new_homomers+old_dusp+old_homomers)
     plt.bar(XS,old_homomers,bottom=new_homomers)

     plt.show(0)
     
class EvolutionResult(object):
    __slots__ = ('L','S_c','mu','dup_rate','asym','discov_types','composition_types','discov_times','div')
    def __init__(self,L=None,S_c=None,mu=None,dup_rate=None):
        self.L = L
        self.S_c = S_c
        self.mu = mu
        self.dup_rate = dup_rate
        self.asym = formTime(L,S_c)/formTime(L//2,S_c)
         
    def addData(self,discov_types, composition_types, discov_times,div):
        self.discov_types = discov_types
        self.composition_types = composition_types
        self.discov_times = discov_times
        self.div = div
     
    def __repr__(self):
        return f'Evolution result for \'L:{self.L}, S_c:{self.S_c}, Dup: {self.dup_rate}\''

def getGammas():
     return {60:1.71, 80:3.32, 100: 3.25, 120:4.46, 140:3.74}

def getRes():
     loadD = (([.83],[0,.05],.00417,60),([.75],[0,.05],.003125,80), ([.74],[0,0.05],.0025,100), ([.7],[0,.05],.0021,120),([.714],[0,.05],.001786,140))
     res=[]
     for d in loadD:
          er, err, ed, (QQ)= loadManyRecords(*d)
          res.extend(QQ)
     return res



          
        
def plotComplexityEvolution(data,renormalise_t0=False,interaction_sets=(10,),N_simulatins=None):

    def deWeight(times,init=1):
        heights = list(np.linspace(1/len(times),init,len(times)))
        unique_times = []
        for index,time in enumerate(reversed(times)):
            if time not in unique_times:
                unique_times.append(time)
            else:
                heights.pop(len(times)-index-1)
        return unique_times[::-1],np.array(heights)

    f, (ax,ax2) = plt.subplots(2)

    for index, c in zip(range(0,len(data),2),('g','m','r','b','c')):
        spline_data = []
        time_window = None

        for sub_I in (1,0):
            raw_times = np.array(sum((data[index+sub_I].discov_times[slice_] for slice_ in interaction_sets),[]))
            time_points, complexity_avg = deWeight(sorted(raw_times),N_simulations or len(interaction_sets))            
            ax.plot(time_points, complexity_avg, ls='--' if sub_I else '-', c=c, label=f'{data[index+sub_I].L},{data[index+sub_I].dup_rate}',alpha=1)

            ## generate spline information for complexity gap
            complexity_spline = InterpolatedUnivariateSpline(time_points, complexity_avg)
            spline_data.append(complexity_spline)

            ##TODO, what is the ideal slice point here?
            time_window = np.linspace(1,time_points[np.argmax(complexity_avg>=3)],250)
            #,complexity_avg[-1]
        
        new_start=formTime(data[index+sub_I].L//2,data[index+sub_I].S_c)/75
        if renormalise_t0:
            time_window = time_window[time_window>=new_start] - new_start
        else:
            ax2.axvline(new_start,c=c)

        complexity_gap = spline_data[0](time_window)-spline_data[1](time_window)
        max_gap = np.argmax(complexity_gap)
        ax.plot([time_window[max_gap]]*2,[spline_data[0](time_window[max_gap]),spline_data[1](time_window[max_gap])],mfc='none',mec=c,marker='o',ms=12,c=c)

        time_window = time_window[:max_gap+1]
        ax2.plot(time_window,spline_data[0](time_window)-spline_data[1](time_window), c=c)

    ax.set_xscale('log')
    ax2.set_xscale('log')
    ax.legend()
    plt.show(block=False)

def plotCompositionBreakdown():
     hets = np.array([[2752,  813],[2234, 1020], [2771,  627], [2427,  960], [3205,  321]],dtype=np.float64)

     scale = np.sum(hets,axis=1)/3500
     print(scale)
     hets = (hets.T/scale).T/3500
     
     asyms= np.array([32.89107322,  8.17588437, 17.44973057,  8.19035052, 38.25540118])

     for L,vals in zip([60,80,100,120,140],hets):
          for val,sgn in zip(vals,[1,-1]):
               plt.plot([sgn*val,0][::sgn],[L]*2)
          print(L,vals,sum(vals))

     plt.plot([0,0],[50,150])
     plt.show(0)


