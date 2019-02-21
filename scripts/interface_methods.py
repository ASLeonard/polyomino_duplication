import numpy as np
from numpy import linalg as LA
import os
from matplotlib import pyplot as plt
from scipy.stats import binom
from math import ceil
from matplotlib.colors import LogNorm
import pandas
import seaborn as sns

from scipy.optimize import curve_fit
def func(x, b):
     return .5*np.exp(-b * x)+.5

def plotHomology(data,L):
     f,(ax1,ax2)=plt.subplots(1,2,sharex=True,sharey=True)
     ax1.set_ylim((0,1))
     ax1.axhline(.5,c='k',lw=3,ls='--')
     ax2.axhline(.5,c='k',lw=3,ls='--')
     cols=['b','r','g','c']
     overalls=[[],[]]
     #GEN=data[0].shape[0]
     for i,sea in enumerate(data):
          mean=np.mean(sea,axis=1)
          GEN=mean.shape[0]
          overalls[i%2].append(mean)
          for j in range(4):
               pass#(ax1 if i%2 else ax2).plot(range(GEN),1-mean[:,j]/L,c=cols[j],alpha=0.5)

     for i in range(2):
          avg=np.mean(np.array(overalls[i]),axis=0)
          for j in range(4):
               if j==0 and i==0:
                    popt, pcov = curve_fit(func, range(GEN), 1-avg[:,j]/L)
                    f=lomax.fit(1-avg[:,])
                    print(popt)
                    (ax1 if i%2 else ax2).plot(range(GEN),func(np.arange(GEN),*popt),c=cols[j],lw=3,ls=':')
               (ax1 if i%2 else ax2).plot(range(GEN),1-avg[:,j]/L,c=cols[j],lw=3)
     plt.show(block=False)
     
def LEH(Y,runs,pop):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/Bomology_Run{}.BIN'.format(runs),dtype=np.uint8).reshape(-1,pop,4)

def loadData(S,fname='Discovs'):
     return [[int(i) for i in line.split()] for line in open('/rscratch/asl47/Discs/{}{:.6f}.BIN'.format(fname,S))]

def loadBinary(S,fname='Discovs',shape=(-1)):
     return np.fromfile('/rscratch/asl47/Discs/{}{:.6f}.BIN'.format(fname,S),dtype=np.uint8).reshape(shape)

def plotDecors(low,high):
     #plt.figure()
     for S in range(low,high+1,2):
          data=loadBinary(S/32,'Decor',(1+(32-S)//2,-1))
          #plt.plot()
          print(S,np.count_nonzero(data,axis=1))

def plotBonds(data):
     f,ax=plt.subplots()

     
     deads=[]
     qq=0
     counters=0
     xc=0
     for run in data:
          
          dead_points=[np.min(np.where(run[:,i]>=24)) if np.any(run[:,i]>=24) else -1 for i in range(3) ]
          deads.append(dead_points)
          if np.argmax(dead_points)==2:
               if -1 in dead_points and dead_points[2]>=0:
                    xc+=1
                    continue
               counters+=1
               #print("wow",dead_points,qq)

               
          qq+=1
          continue
          
     #return deads
     ax.plot(range(100),[24]*100,'k:')
     print("happened ",counters,xc)
     plt.show(block=False)
     
def plotDiscovery(l_I):
     def SF_sym(S_stars):
          return 1/binom(l_I/2,.5).sf(np.ceil(l_I/2*S_stars)-1)
     def SF_asym(S_stars):
          return 1/binom(l_I,.5).sf(np.ceil(l_I*S_stars)-1)
     
     plt.figure()
     s_hats=np.linspace(0,1,65)[34:55]
     
     data=[]
     for s in s_hats:
          for tds,form in zip(loadData(s),('S','A')):
               for td in tds:
                    if td>0:
                         data.append({'td': td, 'sym': form,'thresh':s})

     df=pandas.DataFrame(data)
     ax = sns.violinplot(x="thresh", y="td", hue="sym",data=df, cut=0,gridpoints=1000,palette="Set2", split=True,scale="count", inner="quartile")
     ax.plot(range(len(s_hats)),SF_asym(s_hats),ls='--',c='coral',marker='D',mec='k',mfc='None')
     ax.plot(range(0,len(s_hats),2),SF_sym(s_hats[::2]),ls='--',c='forestgreen',marker='h',mec='w',mfc='None')
     ax.set_yscale('log')
     _=ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
     ax.set(xlabel=r'$\hat{S}_c$', ylabel=r'$\langle \tau_D \rangle$')
     ax.set_title('Discovery time of symmetric and asymmetric interactions')
     
     ab = plt.axes([.3, .6, .2, .2], facecolor='gainsboro')
     for lz,col in zip((32,64,128),('royalblue','darkmagenta','indianred')):
          sts=np.linspace(0,1,lz)[lz//2:]
          l_I=lz
          ab.plot(sts,SF_asym(sts)/SF_sym(sts),ls='--',marker='o',mfc='None',mec=col,c=col,label=lz,lw=2)
     ab.set_yscale('log')
     ab.legend()
     ab.set_title(r'$Q \langle \tau_D \rangle$')
     
     plt.show(block=False)

def plotDecay():
     f,(ax1,ax2) = plt.subplots(1,2)
     s_hats=np.linspace(0,1,65)[34:65]
     data=[np.zeros((int(ceil(len(s_hats)/2)),int(s_hats[0]*32))),np.zeros((len(s_hats),int(s_hats[0]*64)))]

     for i,s in enumerate(s_hats):
          for j,tds in enumerate(loadData(s,'Decays')):
               ar=np.array(tds,dtype=np.float)
               ar[ar==0]=np.nan
               if int(s*64)%2==1 and j%2==0: #is faulty symmetric
                    continue
               if j>0 and j%4==2: #is faulty symmetric gap
                    continue

               data[j%2][i//(2-j%2),(j//2)//(2-j%2)]=np.nanmean(ar)

     rs1=ax1.imshow(data[1].T,norm=LogNorm())
     plt.colorbar(rs1,ax=ax1)
     rs2=ax2.imshow(data[0].T,norm=LogNorm())
     plt.colorbar(rs2,ax=ax2)

     for i,s in enumerate(s_hats):
          ax1.scatter(i,0 if s==1. else 64*(1-RandomWalk(64,0,0,s,1)),marker='D',color='r')
          if i%2==0:
               ax1.scatter(i,0 if s==1. else 64*(1-RandomWalk(32,0,0,s,1)),marker='X',color='c')
          if s==1.:
               exsp=[1]
          else:
               exsp=expected_steps_fast(__matrix(64,s).T)
          max_g=64-int(s*64)
          for j,g in enumerate(range(max_g+1)):
               ax1.text(i,max_g-j,int(exsp[j]),color='w' if int(exsp[j])<=5 else 'k',fontsize=8,va='center',ha='center')

     for i,s in enumerate(s_hats[::2]):
          ax2.scatter(i,0 if s==1. else 32*(1-RandomWalk(32,0,0,s,1)),marker='D',color='r')
          if s==1.:
               exsp=[1]
          else:
               exsp=expected_steps_fast(__matrix(32,s).T)
          max_g=32-int(s*32)
          for j,g in enumerate(range(max_g+1)):
               ax2.text(i,max_g-j,int(exsp[j]),color='w' if int(exsp[j])<=5 else 'k',fontsize=8,va='center',ha='center')
               
     ax1.set_xticks(np.arange(len(s_hats)))
     ax1.set_xticklabels(s_hats,rotation=45)
     

     ax2.set_xticks(np.arange(len(s_hats[::2])))
     ax2.set_xticklabels(s_hats[::2],rotation=45)

     ax1.set_xlabel(r'$\hat{S}_c$')
     ax1.set_ylabel(r'$\Delta S$')

     plt.suptitle('Decay of asymmetrics and symmetrics')
     plt.show(block=False)
     
     
     
BASE_PATH=''
def setBasePath(path,binary_mode):
     if binary_mode:
          default_file='/{}_Run{}.txt'
     else:
          default_file='/{}_Run{}.txt'
          
     global BASE_PATH
     
     BASE_PATH=path+default_file
         
def setLength(length):
     global interface_length
     interface_length=length
     global interface_type
     interface_type={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}[interface_length]

setLength(64)
             
def BindingStrength(base1,base2):
     return 1-bin(np.bitwise_xor(interface_type(base1),reverseBits(base2))).count('1')/interface_length

def reverseBits(value):
    return ~interface_type(int(('{:0'+str(interface_length)+'b}').format(value)[::-1],2))

class Interactions(object):
     __slots__ = ('bonds','strengths')
     def __init__(self,bonds=[],strengths=[]):
          self.bonds=bonds
          self.strengths=strengths
     def __iter__(self):
          for b,s in zip(self.bonds,self.strengths):
               yield (b,s)
     def __repr__(self):
          return "{} interactions".format(len(self.bonds))
         

def T():
     for i in range(20):
          p=LPB(i,250)
          q=np.sum(p,axis=2)
          print(np.where(q==10)[0])
          
def LSHB(run,pop_size):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/Selections_Run{}.txt'.format(run),dtype=np.uint16).reshape(-1,pop_size)

def LPB(run,pop_size):
     return np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/PIDs_Run{}.txt'.format(run),dtype=np.uint8).reshape(-1,pop_size,2)

def LSB(run,pop_size):
     raw_b=np.fromfile('/scratch/asl47/Data_Runs/Bulk_Data/Strengths_Run{}.BIN'.format(run),dtype=np.uint8)
     wheres=np.where(raw_b== 255)[0]
     res=np.empty(wheres.shape,dtype=object)
     low_slice,ith=0,0
     for high_slice in wheres:
          sub_slice=raw_b[low_slice:high_slice]
          if len(sub_slice)==0:
               res[ith]=Interactions()
          else:
               res[ith]=Interactions(tuple(zip(sub_slice[::3],sub_slice[1::3])),1-sub_slice[2::3]/interface_length)
          low_slice=high_slice+1
          ith+=1


     return res.reshape(-1,pop_size)
          
def LoadSelectionHistory(run,S_star,t,mu,gamma):
     selections=[]
     for line in open(BASE_PATH.format('Selections',S_star,t,mu,gamma,run)):
          converted=[int(i) for i in line.split()]
          selections.append(converted)
     return np.array(selections,np.uint16)

def LoadPIDHistory(run):
     phenotype_IDs=[]
     for line in open('/scratch/asl47/Data_Runs/Bulk_Data/PIDs_Run{}.txt'.format(run)):
          phenotype_IDs.append([list(zip(*(iter([int(i) for i in grouping.split()]),) * 2)) for grouping in line.split(',')[:-1]])
          #return converted
          

     return ObjArray(phenotype_IDs)#np.array(phenotype_IDs,dtype=np.uint8)

def LoadStrengthHistory(run,S_star,t,mu,gamma):
     strengths=[]
     for line in open(BASE_PATH.format('Strengths',S_star,t,mu,gamma,run)):
          row=[]
          for species in line.split(',')[:-1]:
               bond_list=[]
               str_list=[]
               for pairing in species.split('.'):
                    if pairing=='':
                         continue
                    pairing=pairing.split()
                    bond_list.append(tuple(int(i) for i in pairing[:2]))
                    str_list.append(1-int(pairing[2])/interface_length)
               row.append(Interactions(bond_list,str_list))
          strengths.append(row)
                
     return ObjArray(strengths)

def LoadPhenotypeTable(run):
     phenotype_table= sorted([[int(i) for i in line.split()] for line in open('/scratch/asl47/Data_Runs/Bulk_Data/PhenotypeTable_Run{}.txt'.format(run))],key=lambda z: z[0])
     return {tuple(px[:2]): tuple(px[2:]) for px in phenotype_table}

def LoadAll(run,params,cwd=None):
     binary_mode=type(params) is int
     setBasePath(cwd if cwd else os.getcwd(),binary_mode)

     if binary_mode:
          s=LSHB(run,params)
          p=LPB(run,params)
          st=LSB(run,params)
     else:     
          st=LoadStrengthHistory(run,*params)
          p=LoadPIDHistory(run,*params)
          s=LoadSelectionHistory(run,*params)
     setBasePath(cwd if cwd else os.getcwd(),True)
     table=LoadPhenotypeTable(run)
     return (s,p,st,table)

def ObjArray(data):
     shape=(len(data),len(data[0]))
     nparr=np.empty(shape,dtype=object)
     nparr[:]=data
     return nparr 

""" DRIFT SECTION """

def __matrix(I_size,S_star):
     N_states=int(I_size*(1-S_star))+1
     val=np.linspace(0,1,I_size+1)[-N_states:]
     rows=[[0]*(N_states+1),[val[0],0,1-val[0]]+[0]*(N_states-2)]
     for i in range(1,N_states-1):
          rows.append([0]*(i)+[val[i],0,1-val[i]]+[0]*(N_states-2-i))
     rows.append([0]*(N_states-1)+[val[-1],0])
     matrix= np.vstack(rows).T
     return matrix
     

def expected_steps_fundamental(Q):
    I = np.identity(Q.shape[0])
    N = np.linalg.inv(I - Q)
    o = np.ones(Q.shape[0])
    return np.dot(N,o)

def expected_steps_fast(Q):
    I = np.identity(Q.shape[0])
    o = np.ones(Q.shape[0])
    return np.linalg.solve(I-Q, o)[1:]-1


def RandomWalk(I_size=64,n_steps=1000,phi=0.5,S_star=0.6,analytic=False):
     s_hats=np.linspace(0,1,I_size+1)
     N=int(I_size*(1-S_star))+1
     
     if analytic:
          analytic_states=__getSteadyStates(__matrix(I_size,S_star)[1:,1:])[1]
          return sum(s_hats[-N:]*analytic_states) if analytic==1 else analytic_states
     
     states=np.array([1]+[0]*(N-1),dtype=float)
     progressive_states=[sum(s_hats[-N:]*states)]

     for i in range(n_steps):
          states=__updateStates(states,s_hats[-N:],phi)
          progressive_states.append(sum(s_hats[-N:]*states))
     return progressive_states

def __updateStates(states,val,phi=0.5):
     states_updating=states.copy()
     for i in range(states.shape[0]):
          states_updating[i]-=states[i]*phi
          if i!=0:
               states_updating[i]+=states[i-1]*phi*(1-val[i-1])
          if i!=states.shape[0]-1:
               states_updating[i]+=states[i+1]*phi*val[i+1]
     return states_updating/sum(states_updating)     

def __getSteadyStates(matrix):
     eigval,eigvec=LA.eig(matrix)
     va=max(eigval)
     ve=eigvec.T[np.argmax(eigval)]
     return va,ve/sum(ve)

