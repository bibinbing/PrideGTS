function mix1 = init_MoG1(x)
% mix1 = init_MoG1(x)
%
% Initialises a 1-d, 1-component Gaussian mixture model 
% for learning using the Variational Bayes framework.
%
% Called from 'init_MoG'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% x                The data vector
%
%
%
% -----------
% Output
% -----------
%
% MIX1 is a data structure with the following fields:
%
% type             'g'
% m                The number of components (=1)
%
%                  In the field priors:
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
% m_0,tau_0        Normal parameters for means
% eta_0            If HMM, Dirichlet parameters for trans.
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
% pi               Mixing coefficients (=1)
% centres          Means
% precs            Precisions
% gammas           Component responsibilities (=1)
%
%
% -------------------------------------------------------
%
% Original code by Will Penny
% Paper: Variational Bayes for 1-dimensional Mixture Models
%        Technical Report (www.robots.ox.ac.uk/~parg)
%
% Modified by Rizwan Choudrey for use in vbICA model
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)


global CHECK_PROGRESS

print_estimates=0;

N=length(x);

mx=mean(x);
vx=var(x);
bx=1/vx;

% PRIORS
% Set component precision priors
mean_precision=bx;
var_precision=bx/N;
b_0=var_precision/mean_precision;
c_0=(mean_precision^2)/var_precision;
% Set component mean priors
m_0=mx;
v_0=(0.3*(max(x)-min(x)))^2; %SHOULD BE OF STANDARD ERROR FORM
tau_0=1/v_0;

% POSTERIORS
b=b_0;
c=c_0;
mm=m_0;
tau=tau_0;
lambda=N;

if CHECK_PROGRESS
  plotMoG(mm,1./(b.*c),lambda./sum(lambda))
  drawnow
end




% Put variables into data structure
mix1.type = 'g';
mix1.m=1;

mix1.priors.lambda_0=N;
mix1.priors.m_0=m_0;
mix1.priors.tau_0=tau_0;
mix1.priors.b_0=b_0;
mix1.priors.c_0=c_0;

mix1.posts.lambda=N;
mix1.posts.mm=mm;
mix1.posts.tau=tau;
mix1.posts.b=b;
mix1.posts.c=c;

mix1.pi=1;
mix1.centres=mm;
mix1.variances=1/(b*c);

mix1.pi=lambda./sum(lambda);
mix1.centres=mm;
mix1.variances=1./(b.*c);
mix1.gammas = ones(size(x));
