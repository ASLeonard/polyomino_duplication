import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap, Normalize


import numpy as np
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
