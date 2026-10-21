function mix1 = init_MoG(x,m,init_method,priors)
% mix1 = init_MoG(x,m,init_method,priors)
%
% Initialises a 1-dimensional Gaussian mixture model
% for learning using the Variational Bayes framework.
%
% Called from 'intialise_mix1d'.
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
%
%
% Optional parameters
%
% init_method   'kmeans','co_mean','rand' (Default = 'kmeans')
% priors        User defined priors
%
%
%
% -----------
% Output
% -----------
%
% MIX1 is a data structure with the following field:
%
% type             'g'
% m                The number of components
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
% pi               Mixing coefficients 
% centres          Means
% precs            Precisions
% gammas           Component posterior probabilities
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


global CHECK_PROGRESS;
x=x(:)';
N=length(x);


if (m==1)
  mix1=init_MoG1(x);
  return;
end

if nargin<3
  init_method='kmeans';
end

if nargin<4
  set_priors = 0;
else
  set_priors = priors.setSource;
end



% PRIORS
% Set mixing priors
lambda_0=5;
if set_priors
  lambda_0 = priors.lambda_0;
end

% Set component mean priors
m_0=mean(x);
v_0=(0.3*(max(x)-min(x)))^2; %SHOULD BE OF STANDARD ERROR FORM
tau_0=1/v_0;

switch init_method
 
 case 'kmeans',

  % Run kmeans on the data
  [kmeans] = kmeans1 (m,x);
  
  %extract component responsibilities
  gammas = kmeans.gammas;

  % Set component precision priors
  kv=kmeans.v;
  mean_precision=mean(1./kv);
  var_precision=std(1./kv)^2;
  b_0=mean(var_precision/mean_precision);
  c_0=mean((mean_precision^2)/var_precision);
  
  % Posterior for means
  lambda=kmeans.pi;
  mm=kmeans.m;
  v=mean(kv)*[1:1:m]./(N/m); % standard error
    
  % Posterior for precisions
  for s=1:m,
    % Set so that b*c=precision and  b^2*c=var_precision
    if (kmeans.v(s)>0)
      % ie. all variances from k-means are non-zero
      precision=(1/kmeans.v(s));
    else
      precision=mean_precision;
    end
    b(s)=var_precision/precision;
    c(s)=(precision^2)/var_precision;
  end
 
 
 case 'co_mean',
  % Cluster on absolute difference from mean
  z=abs(x-mean(x));
  zmix=kmeans1(m,z);

  %extract component responsibilities
  gammas = zmix.gammas;
  
  % Set component precision priors
  zmean=zmix.m.^2;
  mean_precision=1/mean(zmean);
  var_precision=1/std(zmean)^2;
  b_0=mean(var_precision/mean_precision);
  c_0=mean((mean_precision^2)/var_precision);

  % Posterior for means
  mm=mean(x)*ones(1,m);
  v=1./tau_0*ones(1,m);
  
  % Posterior for mixers
  lambda=100*zmix.pi;
  
  % Posterior for precisions
  for s=1:m,
    % Set so that b*c=precision and  b^2*c=var_precision
    precision=1/zmean(s);
    b(s)=var_precision/precision;
    c(s)=(precision^2)/var_precision;
  end
 otherwise
    b_0 = rand*1000;
    c_0 = rand;
    b = repmat(b_0,1,m);
    c = repmat(c_0,1,m);
    gammas = rand(m,N);
    gammsum = repmat(sum(gammas),m,1);
    gammas = gammas./gammsum;
    lambda=(100/m)*ones(1,m);
    rx=max(x)-min(x);
    mm = randn(1,m)*sqrt(rx);
    v = 1./tau_0*ones(1,m);
end


%set user requested priors if necessary
if set_priors
  m_0=priors.m_0;
  b_0=priors.b_0;
  c_0=priors.c_0;
  tau_0=priors.tau_0;
  lambda_0=priors.lambda_0;
end


if CHECK_PROGRESS
  plotMoG(mm,1./(b.*c),lambda./sum(lambda))
  drawnow
end


% Put variables into data structure
mix1.type = 'g';
mix1.m=m;

mix1.priors.lambda_0=lambda_0;
mix1.priors.m_0=m_0;
mix1.priors.tau_0=tau_0;
mix1.priors.b_0=b_0;
mix1.priors.c_0=c_0;

mix1.posts.lambda=lambda;
mix1.posts.mm=mm;
mix1.posts.tau=1./v;
mix1.posts.b=b;
mix1.posts.c=c;


mix1.pi=lambda./sum(lambda);
mix1.centres=mm;
mix1.precs=(b.*c);
mix1.gammas = gammas;


