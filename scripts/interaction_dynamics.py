from math import ceil,floor
from collections import Counter

import matplotlib.pyplot as plt
from scipy.stats import binom, rv_discrete, sem
import numpy as np
from numpy import linalg as LA

##Markov transition matrices
def makeTransitionMatrix(L,S_c):
    N_states = floor(L*round(1-S_c,2))+1
    weights = np.linspace(0,1,L+1)[-N_states:]
    rows = [[0]*(N_states+1), [weights[0],0,1-weights[0]]+[0]*(N_states-2)]
    for i in range(1,N_states-1):
        rows.append([0]*(i)+[weights[i],0,1-weights[i]]+[0]*(N_states-2-i))
    rows.append([0]*(N_states-1)+[weights[-1],0])
    return np.vstack(rows).T

def formationMatrix(L,S_c,drop=0):
    cut = ceil(round(S_c*L,2))+1
    if drop:
        cut = floor(round(S_c*L,2))+2
    m = makeTransitionMatrix(L,0)[1:cut+1,1:cut+1]
    m[-2,-1] = 0
    return m

def dropMatrix(L,S_c):
    fm = formationMatrix(L,1-S_c,1)
    P_m = np.fliplr(np.identity(fm.shape[0]))
    return P_m.dot(fm).dot(P_m)

##solvers
def getSteadyStates(matrix):
    eigval,eigvec = LA.eig(matrix)
    max_val = max(eigval)
    max_vector = eigvec.T[np.argmax(eigval)]
    return max_val,max_vector/sum(max_vector)

def getFundementalMatrix(Q):
    I = np.identity(Q.shape[0])
    return np.linalg.inv(I - Q)

def expectedStepsFundamental(Q):
    return np.dot(getFundementalMatrix(Q),np.ones(Q.shape[0]))

def stepVariance(Q):
    N = getFundementalMatrix(Q)
    t = expectedStepsFundamental(Q)

    return np.dot(2*N - np.identity(Q.shape[0]), t) - t**2
    
##solution
def steadyStates(L,S_c):
    return getSteadyStates(makeTransitionMatrix(L,S_c)[1:,1:])[1]

def formTime(L,S_c):
    bb = binom(L,.5).pmf(range(ceil(round(L*S_c,2))))
    fm = formationMatrix(L,S_c)

    return np.inner(bb,expectedStepsFundamental(fm.T)[:-1]-1)+1

def dropTime(L,S_c):
    dm = dropMatrix(L,S_c)
    times = expectedStepsFundamental(dm.T)[1:]-1
    distrs = getSteadyStates(makeTransitionMatrix(L,S_c)[1:,1:])[1]
    return times.dot(distrs)

def MutualExclusion(n,S_c,L_I=64):
    return (binom(L_I/2,.5).cdf(int(ceil(S_c*L_I/2))-1)**n)*(binom(L_I,.5).cdf(int(ceil(S_c*L_I))-1)**(n*(n-1)/2.))

##approx
def scaler(L,S):
    return np.exp((.7/np.sqrt(L)+.08)*L**(1.35*S))


def loadBinary(S_c,fname,shape=(-1),dtype=np.uint32):
    return np.fromfile(f'{fname}_{S_c:.6f}.BIN',dtype=dtype).reshape(shape)

def loadFormations(S_c,fpath=''):
    return loadBinary(S_c,f'{fpath}Formation',(2,-1))

def loadDecays(L,S_c,fpath=''):
    max_g = ceil(L*(1-S_c)) + 1
    return loadBinary(S_c,f'{fpath}Decay',(max_g,2,-1))

def calculateMeanMetric(data,shape):
    data = data.astype(np.float)
    data[data==0] = np.nan
    return np.nanmean(data,axis=shape)

def loadGammas(L,S_c,fpath='',):
    max_g = ceil(L//2*(1-S_c)) + 1
    return loadBinary(S_c,f'{fpath}Gamma',(max_g,-1),np.int8)

def calculateGammaFactors(L,S_c,fpath=''):
    steady_states = getSteadyStates(makeTransitionMatrix(L//2,S_c)[1:,1:])[1][::-1]
    
    loaded_data = loadGammas(L,S_c,fpath)
    loaded_data = loaded_data.astype(np.float)
    loaded_data[loaded_data==-1] = np.nan
    means = np.nanmean(loaded_data,axis=1)
    sems = sem(loaded_data,axis=1,nan_policy='omit')

    print('Survival percents to first order')
    for g, (m, s) in enumerate(zip(means,sems)):
        print(f'starting strength {(L-g)/L}: {m:.3f}% ±{s:.3f}')
    
    return steady_states.dot(means)*100

def makeDropDistribution(loaded_data,L,S_c):
    distrs = getSteadyStates(makeTransitionMatrix(L,S_c)[1:,1:])[1]

    Target_N = 100000
     
    EMP_set = []

    for index, state in enumerate(distrs[::-1]):
        EMP_set.extend(np.random.choice(loaded_data[index],size=int(state*Target_N)))

    print(np.mean(EMP_set))
    return EMP_set


def makeDiscrete(EMP):
    Cnt = Counter(EMP)

    keys = sorted(list(Cnt.keys()))
    vals = np.array([Cnt[k] for k in keys])
    vals = vals/np.sum(vals)

    cus= rv_discrete(values=(keys,vals))

    return cus

def plotEMP(L,S_c,normed=False,ax=None):
    loaded_data = loadBinary(S_c,'Decay',(-1,100000))

    sym = makeDropDistribution(loaded_data[::4],L//2,S_c)
    asym = makeDropDistribution(loaded_data[1::2],L,S_c)

    if True:
        sym = np.array(sym)
        asym = np.array(asym)

    sym_PD = makeDiscrete(sym)
    asym_PD = makeDiscrete(asym)

    Xs = np.arange(1,101)
    sN = dropTime(L//2,S_c) if normed else 1
    aN = dropTime(L,S_c) if normed else 1
    plt.plot(Xs/sN,sym_PD.pmf(Xs),ls='-.')
    plt.plot(Xs/aN,asym_PD.pmf(Xs),ls='-')

    print(sum(np.all(sym_PD.rvs(size=2) < asym_PD.rvs(size=1)) for _ in range(100000))/100000)

    plt.show(block=False)

dd = [(60,.83),(80,.75),(100,.74),(120,.7),(140,.714)]

def plotExclusion(S_c,Ls,col='orangered',**kwargs):
    xs=np.linspace(1,500,500)
    for L in Ls:
        mut=MutualExclusion(xs,S_c,L)
        plt.plot(xs,mut,c=col,ls='-',**kwargs)
    #plt.plot(xs[:-1],-np.diff(mut),c='royalblue')
    #print -np.diff(mut),sum(-np.diff(mut))
    plt.show(block=False)

def plotInterfaceProbability(l_I,l_g,Nsamps=False):

    def SF_sym(S_stars):
        return binom(l_I/2,.5).sf(np.ceil(l_I/2*S_stars)-1)#*(1./(l_g+1))
    def SF_asym(S_stars):
        return binom(l_I,.5).sf(np.ceil(l_I*S_stars)-1)#-sym(S_stars))/2*((l_g-1.)/(l_g+1))

    def sym_factor(A):
        return float(2)/(A+1)
    def asym_factor(A):
        return float(A-1)/(A+1)

    s_hats=np.linspace(0,1,l_I+1)

    _, ax1 = plt.subplots()
    ax1.plot(s_hats[::2],np.log10(sym_factor(l_g)*SF_sym(s_hats[::2])),ls='',marker='^',c='royalblue')
    ax1.plot(s_hats,np.log10(asym_factor(l_g)*SF_asym(s_hats)),ls='',marker='o',c='firebrick')

    ax2 = ax1.twinx()
     
    ratios=np.log10((sym_factor(l_g)*SF_sym(s_hats))/(asym_factor(l_g)*SF_asym(s_hats)))
    ax2.plot(s_hats,ratios,c='darkseagreen',ls='',marker='h')
    #crossover=np.where(ratios>0)[0][0]
    #ax2.axvline(s_hats[crossover],color='k',ls='--')
    #ax2.axhline(color='k',ls='-',lw=0.2)
    
    Is={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}
    if Nsamps:
        set_length(l_I)
        s_m=np.zeros(l_I+1)
        a_m=np.zeros(l_I+1)
        for _ in range(Nsamps):
            indices=choice(list(cwr(range(l_g),2)))
            if indices[0]!=indices[1]:
                bases=np.random.randint(0,np.iinfo(Is[l_I]).max,dtype=Is[l_I],size=2)
                    
                a_m[np.where(BindingStrength(*bases)>=s_hats)]+=1
            else:
                base=np.random.randint(0,np.iinfo(Is[l_I]).max,dtype=Is[l_I])
                s_m[np.where(BindingStrength(base,base)>=s_hats)]+=1
        s_m2=np.ma.log10(s_m/Nsamps)
        a_m2=np.ma.log10(a_m/Nsamps)
        ax1.plot(s_hats[::2],s_m2[::2],ls='--',c='royalblue')
        ax1.plot(s_hats,a_m2,ls='--',c='firebrick')

    scale_factor=np.log10(asym_factor(l_g)*SF_asym(s_hats))[0]-np.log10(asym_factor(l_g)*SF_asym(s_hats))[-1]
    ax1.text(.2,np.log10(sym_factor(l_g)*SF_sym(.2))-scale_factor*0.03,'symmetric',va='top')
    ax1.text(.2,np.log10(asym_factor(l_g)*SF_asym(.2)+scale_factor*0.05),'asymmetric',va='bottom')

    ax2.text(.1,(ratios[-1]-ratios[0])*.015+ratios[0],'ratio',ha='center',va='bottom')

    ax1.set_ylabel(r'$  \log  Pr $')
    ax2.set_ylabel(r'$\log \mathrm{ratio}$')
    ax1.set_xlabel(r'$\hat{S}$')

    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.show(block=False)
