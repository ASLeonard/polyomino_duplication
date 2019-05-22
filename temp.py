import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
 
def plotPopulation(pids_raw):
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


def add_selection_layer(ax,selections,run,colorize='M'):
    ##add lines in first row
    lines=[[(-.5,i+.5,),(.5,i+.5)] for i in range(selections.shape[1])]

    ##iterate over following generations and add lines
    for g_ind, p_ind in np.ndindex(selections.shape):
        if g_ind==selections.shape[0]:
            continue
        lines.append([(g_ind+.5,selections[g_ind,p_ind]+.5),(g_ind+1.5,p_ind+.5)])

    ##load colour data
    sizes=np.loadtxt('/scratch/asl47/Data_Runs/Bulk_Data/Size_Run{}.txt'.format(run),dtype=np.uint8).reshape(-1)

    
    lc = LineCollection(lines, linewidths=1,linestyle='-',cmap='inferno')
        
    if colorize=='M':
        lc.set_array(mutations)
    elif colorize=='H':
        lc.set_array(homologies)
    elif colorize=='S':
        lc.set_array(sizes)     
            
    plt.colorbar(lc)
    ax.add_collection(lc)

    
from itertools import chain     
def plotPhenotypeFrequencies(pids_raw,run,thresh=0.25):
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
    plt.show(block=False)
