% GNMF_VB_DEMO01		Multiple template, basic VB demo
%
% Illustrates the basic variational Bayes procedure for learning a Poisson model

% Uses :

% Change History :
% Date		Time		Prog	Note
% 21-Jan-2008	 1:28 PM	ATC	Created under MATLAB 6.5.0 (R13)

% ATC = Ali Taylan Cemgil,
% SPCL - Signal Processing and Communications Lab., University of Cambridge, Department of Engineering
% e-mail : atc27@cam.ac.uk


clf

% Number of Rows
W = 40;

% Number of Columns
K = 5;

% Number of templates
I = 3;

% Set prior parameters 
a_tm = 10*ones(W, I);   % Shape
b_tm = ones(W, I);        % Scale
a_ve = ones(I, K);
b_ve = 100*ones(I, K);

% Generate a random template and excitation
% gamrnd is part of statistics toolbox
T = gamrnd(a_tm, b_tm);
V = gamrnd(a_ve, b_ve);

% Generate data
% poissrnd is part of statistics toolbox
x = poissrnd(T*V);

% Execute the fixed point iterations
[g] = gnmf_vb_poisson_mult_fast(x, a_tm, b_tm, a_ve, b_ve, 'EPOCH', 2000, 'UPDATE', 10,...
                                    'tie_a_ve', 'tie_all', ...
                                    'tie_b_ve', 'tie_all',...
                                    'tie_a_tm', 'tie_all', ...
                                    'tie_b_tm', 'tie_all' ...
);


% Plot results
subplot(121)
marginal_plot(x, T, V, 0.8)
%imagesc(x)
title('Data: X=E*V')

subplot(122)
marginal_plot(g.E_T*g.E_V, g.E_T, g.E_V, 0.8)
%imagesc(g.E_T*g.E_V)
title(['E[T]*E[V]   Bound = ' num2str(g.Bound)])


