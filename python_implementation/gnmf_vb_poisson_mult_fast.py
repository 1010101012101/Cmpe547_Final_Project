import numpy as np
import scipy as sp
from scipy import special
from __future__ import division
import numpy.matlib as M

def gnmf_vb_poisson_mult_fast(x,
                            a_tm,
                            b_tm,
                            a_ve,
                            b_ve,
                            EPOCH =1000,
                            METHOD = 'vb',
                            UPDATE = np.inf,
                            t_init = np.random.gamma(a_tm, b_tm/a_tm),
                            v_init = np.random.gamma(a_ve, b_ve/a_ve),
                            tie_a_ve = 'clamp',
                            tie_b_ve = 'clamp',
                            tie_a_tm = 'clamp',
                            tie_b_tm = 'clamp',
                            print_period = 500
                            ):
    W = x.shape[0]
    K = x.shape[1]
    I = b_tm.shape[1]

    M = ~np.isnan(x)
    X = np.zeros(x.shape)
    X[M] = x[M]

    L_t = t_init
    L_v = v_init
    E_t = t_init
    Sig_t = t_init
    Sig_v = v_init

    B = np.zeros(1,EPOCH)
    gammalnX = special.gammaln(X+1)

    for e in range(EPOCH):
        
        LtLv = L_t.dot(L_v)
        tmp = X / (LtLv)
        #check Tranpose
        Sig_t = L_t.dot(tmp.dot(L_v.T))
        Sig_v = L_v.dot(L_t.T.dot(tmp))

        alpha_tm = a_tm + Sig_t
        beta_tm = 1/((a_tm/b_tm) + M.dot(E_v.T))
        E_t = alpha_tm.dot(beta_tm)

        alpha_ve = a_ve + Sig_v
        beta_ve = 1/((a_ve/b_ve) + E_t.T.dot(M))

        E_v = alpha_ve.dot(beta_ve)
