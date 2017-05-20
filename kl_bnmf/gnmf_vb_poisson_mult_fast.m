function [g] = gnmf_vb_poisson_mult_fast(x, a_tm, b_tm, a_ve, b_ve, varargin)
% GNMF_VB_POISSON_MULT_FAST		Variational Bayes for the multiple template model, compact implementation
%
%  [g] = gnmf_vb_poisson_mult(x, a_tm, b_tm, a_ve, b_ve, <parameter, value>)
%
% Inputs :
%	x : Data
%          a_tm, b_tm, a_ve, b_ve : Hyperparameters
%          <parameter,value> : tie_a_ve, tie_b_ve tie_a_tm tie_b_tm : 
%                      'free', : Learn parameters for each cell (prone to overfitting)
%                      'rows', : learn a single parameter for each row 
%                      'cols',: learn a single parameter for each column
%                      'tie_all': learn a single parameter 
%                      'clamp': Don't learn the hyperparameter
%
% Outputs :
%	g : Parameter structure with sufficient statistics
%
% Usage Example : [] = gnmf_vb_poisson_mult_fast();
%
%
% Note	:
% See also

% Uses :

% Change History :
% Date		Time		Prog	Note
% 23-Feb-2008	11:51 AM	ATC	Created under MATLAB 6.5.0 (R13)

% ATC = Ali Taylan Cemgil,
% SPCL - Signal Processing and Communications Lab., University of Cambridge, Department of Engineering
% e-mail : atc27@cam.ac.uk

%% Memory Layout
%% T = T(freq, idx)
%% V = V(idx, time)
  
    
[EPOCH Method Update t_init v_init tie_a_ve tie_b_ve tie_a_tm tie_b_tm print_period] = parse_optionlist(...
    {'EPOCH', 1000;...  % 
     'METHOD', 'vb'; % 
     'UPDATE', Inf; % Update parameters after this epoch
     't_init', gamrnd(a_tm, b_tm./a_tm); 
     'v_init', gamrnd(a_ve, b_ve./a_ve);
     'tie_a_ve', 'clamp'; % {'free', 'rows', 'cols','tie_all', 'clamp'}
     'tie_b_ve', 'clamp'; % {'free', 'rows', 'cols','tie_all', 'clamp'}
     'tie_a_tm', 'clamp'; % {'free', 'rows', 'cols','tie_all', 'clamp'}
     'tie_b_tm', 'clamp'; % {'free', 'rows', 'cols','tie_all', 'clamp'}
     'print_period', 500;
    }, varargin{:});


% Number of frequency bin
W = size(x,1);

% Number of time slices
K = size(x,2);
I = size(b_tm,2);

M = ~isnan(x);
X = zeros(size(x));
X(M) = x(M);

L_t = t_init;
L_v = v_init;
E_t = t_init;
E_v = v_init;
Sig_t = t_init;
Sig_v = v_init;

B = zeros(1, EPOCH);

gammalnX = gammaln(X + 1);

for e=1:EPOCH,
    
    LtLv = L_t*L_v;
    tmp = X./(LtLv);
    Sig_t = L_t.*(tmp*L_v');
    Sig_v = L_v.*(L_t'*tmp);

    alpha_tm = a_tm + Sig_t;
    beta_tm = 1./(a_tm./b_tm + M*E_v' );
    E_t = alpha_tm.*beta_tm;
    
    alpha_ve = a_ve + Sig_v;
    beta_ve = 1./(a_ve./b_ve + E_t'*M);
    E_v = alpha_ve.*beta_ve;

%% --------------------------------------------------------------------------------------------
% Compute the Bound
%if e==EPOCH,
if rem(e, 10)==1, fprintf(1, '*'); end;
if rem(e, print_period)==1 | e==EPOCH,
        g.E_T = E_t;
        g.E_logT = log(L_t);
        g.E_V = E_v;
        g.E_logV = log(L_v);

        g.Bound = -sum(sum(M.*(g.E_T*g.E_V) +  gammalnX )) ...
          + sum(sum( -X.*( ((L_t.*g.E_logT)*L_v + L_t*(L_v.*g.E_logV))./(LtLv) -  log(LtLv) )  )) ...
          + sum(sum(-a_tm./b_tm.*g.E_T - gammaln(a_tm) + a_tm.*log(a_tm./b_tm)  )) ...
          + sum(sum(-a_ve./b_ve.*g.E_V - gammaln(a_ve) + a_ve.*log(a_ve./b_ve)  )) ...
          + sum(sum( gammaln(alpha_tm) + alpha_tm.*(log(beta_tm) + 1)  )) ...
          + sum(sum( gammaln(alpha_ve)  + alpha_ve.*(log(beta_ve) + 1)  ));
        
        g.a_ve = a_ve;
        g.b_ve = b_ve;
        g.a_tm = a_tm;
        g.b_tm = b_tm;
        
        fprintf(1, '\nBound = %f\t a_ve = %f \t b_ve = %f \t a_tm = %f \t b_tm = %f\n', g.Bound, a_ve(1), b_ve(1), a_tm(1), b_tm(1));

    end
    %% --------------------------------------------------------------------------------------------
    if e==EPOCH,
        break
    end;
    
    L_t = exp(psi(alpha_tm)).*beta_tm;
    L_v = exp(psi(alpha_ve)).*beta_ve;

    if e>Update,
    if ~strcmp( tie_a_tm, 'clamp'),
        Z = E_t./b_tm - (log(L_t) - log(b_tm));
        switch tie_a_tm,
          case 'free',
            a_tm = gnmf_solvebynewton(Z, a_tm);
          case 'rows',
            a_tm = gnmf_solvebynewton(sum(Z,1)/W, a_tm);
          case 'cols',
            a_tm = gnmf_solvebynewton(sum(Z,2)/I, a_tm);
          case 'tie_all',
            a_tm = gnmf_solvebynewton(sum(Z(:))/(W*I), a_tm);
            % case 'clamp', do nothing
        end;
    end;

    switch tie_b_tm,
      case 'free',
        b_tm = E_t;
      case 'rows',
        b_tm = repmat(sum(a_tm.*E_t, 1)./sum(a_tm,1), [W 1]);
      case 'cols',
        b_tm = repmat(sum(a_tm.*E_t, 2)./sum(a_tm,2), [1 I]);
      case 'tie_all',
        b_tm = sum(sum(a_tm.*E_t))./sum(a_tm(:)).*ones(W, I);
        % case 'clamp', do nothing
    end;

    if ~strcmp( tie_a_ve, 'clamp'),
            Z = E_v./b_ve - (log(L_v) - log(b_ve));
            switch tie_a_ve,
              case 'free',
                a_ve = gnmf_solvebynewton(Z, a_ve);
              case 'rows',
                a_ve = gnmf_solvebynewton(sum(Z,1)/I, a_ve);
              case 'cols',
                a_ve = gnmf_solvebynewton(sum(Z,2)/K, a_ve);
              case 'tie_all',
                a_ve = gnmf_solvebynewton(sum(Z(:))/(I*K), a_ve);
                % case 'clamp', do nothing
            end;
        end;

    switch tie_b_ve,
      case 'free',
        b_ve = E_v;
      case 'rows',
        b_ve = repmat(sum(a_ve.*E_v, 1)./sum(a_ve,1), [I 1]);
      case 'cols',
        b_ve = repmat(sum(a_ve.*E_v, 2)./sum(a_ve,2), [1 K]);
      case 'tie_all',
        b_ve = sum(sum(a_ve.*E_v))./sum(a_ve(:)).*ones(I, K);
        % case 'clamp', do nothing
    end;


    
    end;

    %    if rem(e, 100)==1,
    %    fprintf(1, 'a_ve = %f\n', a_ve(1));
    %end;
end;

