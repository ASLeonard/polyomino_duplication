import numpy as np
from numpy import linalg as LA
from math import ceil,floor
from scipy.stats import binom

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
    s_hats = np.linspace(0,1,L+1)
    N = floor(L*(1-S_c))+1
     
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

##approx
def scaler(L,S):
     return np.exp((.7/np.sqrt(L)+.08)*L**(1.35*S))
