function network  = initialise_ica(data,source_type,states,isonoise,eta)
% network  = initialise_ica(data,src_type,states,isonoise)
% 
% Initialisation of vbICA model parameters. Called from 'initialise_ica'.
% 
%
% -----------
% Input
% -----------
%
% Necessary paramteres
%
% data        No. rows = No. sensors, No. columns = No. data points
% source_type Mixture of 'g'aussians     (Default)
%             Mixture of 'e'xponentials
%             Mixture of 't'runcated gaussians
%             'hg','he','ht' - HMM versions
% states      No. states per source model eg [2 4 3 5] for 4 sources
%
%
% Optional parameters
%
% isonoise    1: learn isotropic noise  (Default)
%             0: learn non-isotropic diagonal noise
% eta         Class responsibilities if being trained as part of a 
%             Mixture of ICAs, otherwise = 1
%
%
%
% -----------
% Output
% -----------
%
% NETWORK is a data structure with the following fields:
%
% source_type   The source type
% source_states The source state vector
% sources       The cell array of source models
% alphas        The initialised ARD parameters
% noise         The initialised noise precision (precision of e).
% mean          The initialised noise mean
% priors        The model priors
% posts         The model posteriors
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)


if nargin<4
  isonoise = 1;
end

if nargin<5
  eta=1;
end

HMM=0;
if length(source_type)==2
  src_type = source_type(2);
  HMM=1;
else
  src_type = source_type;
end



if src_type == 'g'
    POS_SRC = 0;
else
    POS_SRC = 1;
end


%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


sources = length(states);
[sensors,points] = size(data);



%---------------------------------set priors------------------------------------
b_alpha_0 = 1000;
c_alpha_0 = 1e-3;
b_Lam_0 = 1000;
c_Lam_0 = 1e-3;
priors.mix.b_alpha_0 = b_alpha_0;
priors.mix.c_alpha_0 = c_alpha_0;
priors.noise.b_Lam_0 = b_Lam_0;
priors.noise.c_Lam_0 = c_Lam_0;

switch src_type
 case 'g'
  priors.noise.mb0 = 1;
  priors.noise.mn0 = 0;
  priors.source.m_0 = 0;
  priors.source.tau_0 = 1;
  priors.source.b_0 = 1000;
  priors.source.c_0 = 1e-3;
  priors.source.lambda_0 = 5;
  priors.source.setSource = 0;
  net_init = 'SVD';
  source_init = 'kmeans';
 case 'e'
  priors.noise.mb0 = inf;
  priors.noise.mn0 = 0;
  priors.source.m_0 = 0;
  priors.source.tau_0 = inf;
  priors.source.b_0 = 1000;
  priors.source.c_0 = 1e-3;
  priors.source.lambda_0 = 5;
  priors.source.setSource = 0;
  net_init = 'rand';
  source_init = 'rand';
 case 't'
  priors.noise.mb0 = inf;
  priors.noise.mn0 = 0;
  priors.source.m_0 = 0;
  priors.source.tau_0 = inf;
  priors.source.b_0 = 1000;
  priors.source.c_0 = 1e-3;
  priors.source.lambda_0 = 5;
  priors.source.setSource = 0;
  net_init = 'rand';
  source_init = 'rand';
end  
if HMM
  priors.source.eta0=5;
  priors.source.eta_ii=5;
end
%-------------------------------set priors-----------------------------------




%-----------------------------set posteriors---------------------------------
[lambda_hat,m_H_hat,mu_x_hat,source,mn] = init_vbica(data,source_type,states,isonoise,priors,net_init,source_init,eta);




mb = repmat(priors.noise.mb0,sensors,1);
beta_x_hat = 1./((1e-4)*ones(size(mu_x_hat)));

%initialise ARD
b_alpha_hat = b_alpha_0*ones(sources,1);
c_alpha_hat = c_alpha_0;
alpha = b_alpha_hat*c_alpha_hat;
alpha_hat = repmat(alpha',sensors,1);
%-----------------------------set posteriors---------------------------------


%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




% Put variables into data structure
network.source_type = source_type;
network.source_states = states;
network.sources = source;
network.alphas = alpha;
network.noise = lambda_hat;
network.mean = mn;

network.priors = priors;
network.posts.x_mean=mu_x_hat;
network.posts.x_prec=beta_x_hat;
network.posts.mix.H_mean=m_H_hat;
network.posts.mix.H_prec=alpha_hat;
network.posts.alpha.b=b_alpha_hat;
network.posts.alpha.c=c_alpha_hat;
if isonoise
    c=points*sensors/2;
else
    c=points/2;
end
network.posts.noise.b = lambda_hat./c;
network.posts.noise.c = c;
network.posts.noise.mn = mn;
network.posts.noise.mb = mb;

