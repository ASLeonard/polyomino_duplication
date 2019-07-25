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

def plotRidge(df):
     sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

     # Create the data
     #rs = np.random.RandomState(1979)
     #x = rs.randn(500)
     #g = np.tile('H', 50)
     #df = pd.DataFrame(dict(x=x, g=g))
     #m = df.g.map(ord)
     #df["x"] += m

     # Initialize the FacetGrid object
     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
     g = sns.FacetGrid(df, row="occurence", hue="occurence", aspect=15, height=.5, palette=pal,sharey=False)
     
     # Draw the densities in a few steps
     g.map(sns.kdeplot, "homology", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
     g.map(sns.kdeplot, "homology", clip_on=False, color="k", lw=2, bw=.2)
     g.map(plt.axhline, y=0, lw=2, clip_on=False)
     g.map(plt.axvline, x=64, lw=2, clip_on=False)
     
     
     # Define and use a simple function to label the plot in axes coordinates
     def label(x, color, label):
          ax = plt.gca()
          ax.text(0, .2, label, fontweight="bold", color=color,
                  ha="left", va="center", transform=ax.transAxes)
          

     g.map(label, "occurence")
          
     # Set the subplots to overlap
     g.fig.subplots_adjust(hspace=-.25)
     
     # Remove axes details that don't play well with overlap
     g.set_titles("")
     g.set(yticks=[])
     g.despine(bottom=True, left=True)
     plt.show(block=False)
