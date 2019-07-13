import numpy as np
from numpy import linalg as LA
from math import ceil,floor
from scipy.stats import binom

##Markov transition matrices
def __matrix(I_size,S_star):
    N_states=int(I_size*(1-S_star))+1
    val=np.linspace(0,1,I_size+1)[-N_states:]
    rows=[[0]*(N_states+1),[val[0],0,1-val[0]]+[0]*(N_states-2)]
    for i in range(1,N_states-1):
        rows.append([0]*(i)+[val[i],0,1-val[i]]+[0]*(N_states-2-i))
    rows.append([0]*(N_states-1)+[val[-1],0])
    matrix= np.vstack(rows).T
    return matrix

def formation_matrix(L,S_c,drop=0):
    cut = ceil(round(S_c*L,2))+1
    if drop:
        cut = floor(round(S_c*L,2))+2
    print(cut)
    m = __matrix(L,0)[1:cut+1,1:cut+1]
    m[-2][-1] = 0
    return m

def drop_matrix(I_size,S_star):
    fm = formation_matrix(I_size,1-S_star,1)
    P_m = np.fliplr(np.identity(fm.shape[0]))
    return P_m.dot(fm).dot(P_m)

##solvers
def __getSteadyStates(matrix):
    eigval,eigvec=LA.eig(matrix)
    va=max(eigval)
    ve=eigvec.T[np.argmax(eigval)]
    return va,ve/sum(ve)


def expected_steps_fundamental(Q):
    I = np.identity(Q.shape[0])
    N = np.linalg.inv(I - Q)
    o = np.ones(Q.shape[0])
    
    return np.dot(N,o)

##solutions
def getExp(L,S_c):
    bb = binom(L,.5).pmf(range(ceil(round(L*S_c,2))))
    fm = formation_matrix(L,S_c)
    return np.inner(bb,expected_steps_fundamental(fm.T)[:-1]-1)+1

def steady_states(L,S_c):
    s_hats=np.linspace(0,1,L+1)
    N=int(L*(1-S_c))+1
     
    analytic_states=__getSteadyStates(__matrix(L,S_c)[1:,1:])[1]
    return analytic_states

def drop_time(L,S_c):
    ##fix rounding error
    dm = drop_matrix(L,S_c)
    times = expected_steps_fundamental(dm)[1:]
    distrs = __getSteadyStates(__matrix(L,S_c)[1:,1:])[1]
    return times.dot(distrs)

##approx
def scaler(L,S):
     return np.exp((.7/np.sqrt(L)+.08)*L**(1.35*S))
