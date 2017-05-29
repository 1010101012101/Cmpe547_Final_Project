import numpy as np
import scipy as sp
from scipy import special
import numpy.matlib as M

W = 40
K = 5
I = 3

a_tm = 10 * np.ones([W,I])
b_tm = np.ones([W,I])
a_ve = np.ones([I,K])
b_ve = 100 * np.ones([I,K])

T = np.random.gamma(a_tm,b_tm)
V = np.random.gamma(a_ve,b_ve)

x = np.random.poisson(T.dot(V))

hoho = gnmf_vb_poisson_mult_fast(x,a_tm,b_tm,a_ve,b_ve,
                                EPOCH=2000,
                                Update =10,
                                tie_a_ve='tie_all',
                                tie_b_ve='tie_all',
                                tie_a_tm='tie_all',
                                tie_b_tm='tie_all')
