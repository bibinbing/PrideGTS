function mix1 = initialise_mix1d(x,m,source_type,init_method,priors)
% mix1 = initialise_mix1d(x,m,source_type,init_method,priors)
%
% Initialises a 1-dimensional mixture model for
% learning using the Variational Bayes framework.
%
% Called by 'mixmodel1d' and 'init_vbica'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% x             The data vector 
% m             Number of mixture components
% source_type   'g', 'e', 't' - Mixture of Gaussians, 
%               exponentials or truncated Gaussians
%               'hg','he','ht' - HMM versions
%
%
% Optional parameters
%
% init_method   'co_mean', 'kmeans', 'rand' for MoG
%               'rand' only for MoP
% priors        User defined priors
%
%
%
% -----------
% Output
% -----------
%
% The fields in MIX1 are:
%
% type             See source_type above
% m                The number of components
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
% gammas           Component posterior probabilities
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




global CHECK_PROGRESS;

if nargin<4
  init_method='rand';
end

if nargin<5
  set_priors = 0;
  PRIORS = 0;
else
  set_priors = priors.setSource;
  PRIORS = 1;
end


x=x(:)';
N=length(x);

if length(source_type)==2
  src_type = source_type(2);
  HMM=1;
else
  src_type = source_type;
  HMM=0;
end


switch src_type
 case 'g'
  if PRIORS
    mix1 = init_MoG(x,m,init_method,priors);
  else
    mix1 = init_MoG(x,m,init_method);
  end
  
 otherwise
  if PRIORS
    mix1 = init_MoP(x,m,src_type,priors);
  else
    mix1 = init_MoP(x,m,src_type);
  end
end



%Extra prior needed for HMM transition matrix
if HMM

  lambda_0 = 0.5;
  eta0 = 5;
  eta_ii = 5;
  
  if set_priors
    lambda_0 = priors.lambda0;
    eta0 = priors.eta0;
    eta_ii = priors.eta_ii;
  end

  eta_0 = repmat(eta0,m,m);
  sumeta = repmat(sum(eta_0,2),1,m);
  mix.priors.lambda_0 = lambda_0;
  mix1.type = source_type;
  mix1.priors.eta_0 = eta_0+eta_ii*eye(m);
  mix1.posts.eta = eta_0;
  mix1.P=eta_0./sumeta;
end





