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
                            Method = 'vb',
                            Update = np.inf,
                            t_init = np.random.gamma(a_tm, b_tm/a_tm),
                            v_init = np.random.gamma(a_ve, b_ve/a_ve),
                            tie_a_ve = 'clamp',
                            tie_b_ve = 'clamp',
                            tie_a_tm = 'clamp',
                            tie_b_tm = 'clamp',
                            print_period = 500
                            ):

    # Result initialiation
    g = dict()
    g['E_T'] = None
    g['E_logT'] = None
    g['E_V'] = None
    g['E_logV'] = None
    g['Bound'] = None
    g['a_ve'] = None
    g['b_ve'] = None
    g['a_tm'] = None
    g['b_tm'] = None

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
        # Compute the bound
        if(e%10 == 10):
            print("*")
        if(e%print_period == 1 or e == EPOCH):
            g['E_T'] = E_t
            g['E_logT'] = np.log(L_t)
            g['E_V'] = E_v
            g['E_logV'] = np.log(L_v)

            g['Bound'] = -np.sum(np.sum(M * (g['E_V'].dot(g['E_V'])) + gammalnX))
                        + np.sum(np.sum(-X * ( ((L_t * g['E_logT']).dot(L_v) + L_t.dot(L_v * g['E_logV']))/(LtLv) - np.log(LtLv) ) ))
                        + np.sum(np.sum((-a_tm/b_tm)* g['E_T'] - special.gammaln(a_tm) + a_tm * np.log(a_tm /b_tm)))
                        + np.sum(np.sum((-a_ve/b_ve)* g['E_V'] - special.gammaln(a_ve) + a_ve * np.log(a_ve /b_ve)))
                        + np.sum(np.sum( special.gammaln(alpha_tm) + alpha_tm * np.log(beta_tm) + 1))
                        + np.sum(np.sum(special.gammaln(alpha_ve) + alpha_ve * np.log(beta_ve) + 1 ))

            g['a_ve'] = a_ve
            g['b_ve'] = b_ve
            g['a_tm'] = a_tm
            g['b_tm'] = b_tm

            print( '\nBound = %f\t a_ve = %f \t b_ve = %f \t a_tm = %f \t b_tm = %f\n', g.Bound, a_ve[0], b_ve[0], a_tm[0], b_tm[0])
        if (e == EPOCH):
            break;
        L_t = np.exp(special.psi(alpha_tm)) * beta_tm
        L_v = np.exp(special.psi(alpha_ve)) * beta_ve

        Z = None
        if( e> Update):
            if(not tie_a_tm == 'clamp' ):
                Z = (E_t / b_tm) - (np.log(L_t) - np.log(b_tm))
                if(tie_a_tm == 'clamp'):
                    a_tm = gnmf_solvebynewton(Z,a_tm)
                elif(tie_a_tm == 'rows'):
                    a_tm = gnmf_solvebynewton(np.sum(Z,0)/W, a_tm)
                elif(tie_a_tm == 'cols'):
                    a_tm = gnmf_solvebynewton(np.sum(Z,1)/I, a_tm)
                elif(tie_a_tm == 'tie_all'):
                    a_tm = gnmf_solvebynewton(np.sum(Z)/(W.dot(I)), a_tm)

            if(tie_b_tm == 'free'):
                b_tm = E_t
            elif(tie_b_tm == 'rows'):
                b_tm = M.repmat(np.sum(a_tm * E_t,0)/np.sum(a_tm,0),W,1)
            elif(tie_b_tm == 'cols'):
                b_tm = M.repmat(np.sum(a_tm * E_t,1)/np.sum(a_tm,1),1,I)
            elif(tie_b_tm == 'tie_all'):
                b_tm = (np.sum(a_tm*E_t)/ np.sum(a_tm)) * np.ones([W,I])

            if(not tie_a_v2 == 'clamp' ):
                Z = (E_v / b_ve) - (np.log(L_v) - np.log(b_ve))
                if(tie_a_ve == 'clamp'):
                    a_ve = gnmf_solvebynewton(Z,a_ve)
                elif(tie_a_ve == 'rows'):
                    a_ve = gnmf_solvebynewton(np.sum(Z,0)/I, a_ve)
                elif(tie_a_ve == 'cols'):
                    a_ve = gnmf_solvebynewton(np.sum(Z,1)/K, a_ve)
                elif(tie_a_ve == 'tie_all'):
                    a_ve = gnmf_solvebynewton(np.sum(Z)/(I.dot(K)), a_ve)

            if(tie_b_ve == 'free'):
                b_ve = E_v
            elif(tie_b_ve == 'rows'):
                b_ve = M.repmat(np.sum(a_ve * E_v,0)/np.sum(a_ve,0),I,1)
            elif(tie_b_tm == 'cols'):
                b_ve = M.repmat(np.sum(a_ve * E_v,1)/np.sum(a_ve,1),1,K)
            elif(tie_b_tm == 'tie_all'):
                b_ve = (np.sum(a_ve*E_v)/ np.sum(a_ve)) * np.ones([I,K])
return g
