from math import ceil,floor
from collections import Counter

import matplotlib.pyplot as plt
from scipy.stats import binom, rv_discrete
import numpy as np
from numpy import linalg as LA

##Markov transition matrices
def makeTransitionMatrix(L,S_c):
    N_states = floor(L*(1-S_c))+1
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


def loadBinary(S,fname='Discovery',shape=(-1)):
    return np.fromfile('{}_{:.6f}.BIN'.format(fname,S),dtype=np.uint32).reshape(shape)
## testing

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

def calcGamma(L,S_c):
    sx= getSteadyStates(makeTransitionMatrix(L//2,S_c)[1:,1:])[1][::-1]
    loaded_data = loadBinary(S_c,'Decay',(-1,100000))
    print('M',np.mean(loaded_data,axis=1))
    print('S',np.std(loaded_data,axis=1))
    return sx.dot(np.mean(loaded_data,axis=1))*100

dd = [(60,.83),(80,.75),(100,.74),(120,.7),(140,.714)]



def plotExclusion(S_c,Ls,col='orangered'):
    xs=np.linspace(1,500,500)
    for L in Ls:
        mut=MutualExclusion(xs,S_c,L)
        plt.plot(xs,mut,c=col,marker='h',ls='')
    #plt.plot(xs[:-1],-np.diff(mut),c='royalblue')
    #print -np.diff(mut),sum(-np.diff(mut))
    plt.show(block=False)
