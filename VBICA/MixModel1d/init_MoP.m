function mix1 = init_MoP(x,m,source_type,priors)
% mix1 = init_MoP(x,m,source_type,priors)
%
% Initialises a 1-dimensional positive mixture model
% for learning using the Variational Bayes framework.
%
% Called from 'initiaise_mix1d'.
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
% source_type   'e', 't' - Mixture of exponentials,
%                or truncated Gaussians
%
%
% Optional parameters
%
% priors        User defined priors
%
%
%
% -----------
% Output
% -----------
%
% MIX1 is a data structure with the following fields:
%
% type             'e' or 't'
% m                The number of components
%
%                  In the field 'priors', parameters for
%                  prior distributions:
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
% m_0,tau_0        Normal parameters for means (dummy)
% eta_0            If HMM, Dirichlet parameters for trans.
%                  matrix
%
%                  In the field 'posts', parameters for
%                  posterior distributions:
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% mm,v             Normal parameters for means (dummy)
% eta              If HMM, Dirichlet parameters for trans.
%                  matrix
%
%                  In the field posts:
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% mm,v             Normal parameters for means
% eta              If HMM, Dirichlet parameters for trans.
%                  matrix
%
%                  Expected posterior values:
% pi               Mixing coefficients 
% centres          Means
% precs            Precisions
% gammas           Component posterior probabilities
%
%
% -------------------------------------------------------
%
% Original code by Rizwan Choudrey
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)


global CHECK_PROGRESS;
x=x(:)';
N=length(x);


if nargin<4
  set_priors = 0;
else
  set_priors = priors.setSource;
end



% PRIORS
lambda_0=5;
m_0 = 0;
b_0 = 1000;
c_0 = 1e-3;
tau_0 = inf;
mm = 0;
tau = inf;
z=[x(1:2:N),-x(2:2:N)];
mean_prec = 1/var(z');

%set user requested priors if necessary
if set_priors
  b_0 = priors.b_0;
  c_0 = priors.c_0;
  tau_0=priors.tau_0;
  lambda_0=priors.lambda_0;
end


lambda = repmat(lambda_0,1,m);
precs = rand(1,m);
gammas = rand(m,N);
gammas = gammas./repmat(sum(gammas,1),m,1);
c = repmat(N./(2*m),1,m);
b = precs./c;

if CHECK_PROGRESS
  switch source_type
   case 't'
    plotRMoG(1./(b.*c),lambda./sum(lambda),max(x)*10)
    drawnow
   case 'e'
    plotMoE(1./(b.*c),lambda./sum(lambda),max(x)*10)
    drawnow
  end
end


% Put variables into data structure
mix1.type = source_type;
mix1.m=m;

mix1.priors.lambda_0=lambda_0;
mix1.priors.m_0=m_0;
mix1.priors.tau_0=tau_0;
mix1.priors.b_0=b_0;
mix1.priors.c_0=c_0;

mix1.posts.lambda=lambda;
mix1.posts.mm=mm;
mix1.posts.tau=tau;
mix1.posts.b=b;
mix1.posts.c=c;

mix1.pi=lambda./sum(lambda);
mix1.centres=mm;
mix1.precs=(b.*c);   
mix1.gammas = gammas;

