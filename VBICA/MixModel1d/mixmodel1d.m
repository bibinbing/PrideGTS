function mix1 = mixmodel1d(data,comps,source_type,tol,max_steps, ...
			   DRAW)
% mix1 = mixmodel1d(data,comps,source_type,tol,max_steps,DRAW)
%
% Train a 1-dimensional mixture model or HMM using the 
% Variational Bayes framework.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% data         The observation vector
% comps        Number of components in mixture
% source_type  'g', 'e', 't' - Mixture of Gaussians, 
%              exponentials or truncated Gaussians
%              'hg','he','ht' - HMM versions
%
%
% Optional parameters
%
% tol       Convergence tolerance          (Default = 1e-5)
% max_steps Max number of iteration steps  (Default = 500)
% DRAW      Plot pdf at every iteration    (Default = 0)
%
%
% Note:     Exponentials take longer to converge, particularly 
%           if source pdfs not well described by mixture of   
%           exponentials. Therefore, tol=1e-6 for MoE. Truncated 
%           Gaussians are more robust and have better convergence.
%
%
%
% -----------
% Output
% -----------
%
% MIX1 is a data structure with the following fields:
%
% type             'g', 'e', 't' - Mixture of Gaussians, 
%                  exponentials or truncated Gaussians
%                  'hg','he','ht' - HMM versions
% m                The number of comp
% f_hidd           The negative free energy of hidden vars.
% f_params         The negative free energy of the params.
% fm               The total negative free energy
%
%                  In the field priors:
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
% m_0,tau_0        Normal parameters for means
%
%                  In the field posts:
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% mm,tau           Normal parameters for means
%
%                  Expected posterior values:
% pi               Mixing coefficients
% centres          Means 
% precs            Precisions
% gammas           Component probabilities
%
%
% If HMM, also
%
% pi               Initial state probabilities
% P                Calculated trans. probs. for HMM
% posts.eta        Dirichlet parameters for P
% posts.lambda     Dirichlet parameters for pi
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

addpath ../utils

global CHECK_PROGRESS;


if nargin<4
  tol = 1e-5;
end

if nargin<5
  max_steps = 500;
end

if nargin<6
  CHECK_PROGRESS = 0;
else
  CHECK_PROGRESS = DRAW;
end


if (length(source_type)==2)
  src_type = source_type(2);
else
  src_type = source_type;
end

switch src_type
 case 'g'
  init_method = 'kmeans';
% init_method = 'co_mean';
 case 't'
  init_method = 'rand';
 case 'e'
  init_method = 'rand';
  tol = 1e-6; % 'e' requires longer convergence
end

mix1 = initialise_mix1d(data,comps,source_type,init_method);
mix1 = learn_mix1d(mix1,data,data.^2,tol,max_steps);