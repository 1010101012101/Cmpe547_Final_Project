
This folder contains a Matlab implementation of variational Bayes for KL (Kullback-Leibler) divergence based Non Negative Matrix Factorisation, as discussed in [1]. 

To run a demo, simply extract all files to a directory, in Matlab change to the directory and run the demo: 

>>  gnmf_vb_demo01

This generates a positive matrix X from the NMF generative model, than computes a positive decomposition with hyperparameter adaptation.

Note: This implementation needs the statistics toolbox of Matlab. Actually, we only use the Poisson and Gamma random variable generators poissrnd and gamrnd. 

File List:

gnmf_vb_poisson_mult_fast.m            NMF implementation with bound computation
gnmf_solvebynewton.m                       Hyper parameter optimisation
marginal_plot.m                                  Visualisation Utility
parse_optionlist.m                               Parameter parsing utility



Ref:
[1] A. T. Cemgil. Bayesian inference in non-negative matrix factorisation models. Technical Report CUED/F-INFENG/TR.609, University of Cambridge, July 2008. Accepted for publication to Computational Intelligence and Neuroscience

%11-May-2009	 3:09 PM	ATC	
