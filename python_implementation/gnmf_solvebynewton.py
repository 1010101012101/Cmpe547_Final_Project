from __future__ import division
import numpy as np
import scipy as sp
from scipy import special
import numpy.matlib as M

def gnmf_solvebynewton(c, a0 = None):

    if a0 is None:
        a0 = 0.1 * np.ones(np.shape(c))

    M, N = np.shape(a0)
    if len(np.shape(c)) == 0:
        Mc , Nc = 1,1
    else:
        Mc, Nc = np.shape(c)



    a = None
    cond = 0

    if (M == Mc and N == Nc):
        a = a0
        cond = 1

    elif (Mc == 1 and Nc >1):
        cond = 2
        a = a0[0,:]
    elif (Mc > 1 and Nc == 1):
        cond = 3
        a = a0[:,0]
    elif (Mc == 1 and Nc == 1):
        cond = 4
        a = a0[0,0]

    a2 = None
    for index in range(10):
        a2 = a - (np.log(a) - special.polygamma(0,a) + 1 - c)/
        (1/a - special.polygamma(1,a))
        idx = np.where(a2<0)[0]
        if( True in idx):
            a2[a2<0] = a2[a2<0] / 2
        a = a2

    if(cond == 2):
        a = M.repmat(a,M,1)
    elif(cond == 3):
        a = M.repmat(a,1,N)
    elif(cond == 4):
        a = a * np.ones([M,N])

    return a
