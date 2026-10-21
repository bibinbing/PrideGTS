function recon = unmix2(data,net,tol,max_steps,POSITIVE,old_x,eta);
% recon = unmix1(data,net,tol,max_steps,POSITIVE);
% 
% Use model trained by vbICA1 to unmix data and infer hidden
% source signals and source states.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        No. rows = No. sensors, No. columns=No. data points
% net         ICA network to use for unmixing
%
%
% Optional Parameters:
%
% tol         Convergence tolerance (default=1e-4)  
% max_steps   Max no_steps to cease (default=500) 
% POSITIVE    Positive-only source signals (default=1 for 'e', 't')
%                                          (default=0 for 'g')
%
%
%
% -----------
% Output
% -----------
%
% RECON is a data structure with the following fields:
%
% signals     Inferred hidden source signals
% means       Means of hidden source signal posterior
% precs       Precisions of hidden source signal posterior
% state_probs Probabilities of hidden source states
% state_paths Most likely state paths (Viterbi if HMM)
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
% function [x,x_prec,fx,fx_d] = unmix2(data,net,tol,max_steps,eta);
% 
% Use model trained by vbICA2 to unmix dataset.
%
%
% Necessary parameters
%
% data        No_rows=No_sensors.
%             No_columns=No_data points.
% net         ICA network to use for unmixing.
%
%
%
% Optional Parameters
%
% tol         Convergence tolerance. (default=1e-4)  
% max_steps   Max no_steps to cease.   
% eta         Class responsibilities if being trained as part of a Mixture of ICAs.
%             Otherwise = 1;
%
%------------------------------------------------------------------------------------
% Output
%
% x           Expectations of hidden source vectors under posterior.
% x_prec      Precision of expectations of hidden source vectors under posterior.
% fx          Energy of reconstruction.
%
% Optional Parameters
%
% fx_d        Energy of reconstruction per data point.

addpath ../utils ../Shared ../MixModel1d  

%--------------------Error catching and default setting--------------------
if nargin<3
  tol = 1e-4;
end

if nargin<4
  max_steps = 500;
end


if length(net.source_type)==2
  HMM = 1;
  src_type = net.source_type(2);
else
  HMM = 0;
  src_type = net.source_type;
end

if nargin<5
  switch src_type
   case 'g'
    POSITIVE = 0;
   case 'e'
    POSITIVE = 1;
   case 't'
    POSITIVE = 1;
  end
end

if nargin<6 | isempty(old_x)
  init_x = 1;
else
  init_x = 0;
end


if nargin<7
  eta = 1;
end

if tol<eps
  error('Please choose positive TOL!')
end

if max_steps<1
  error('Please choose MAX_STEPS > 1!')
end

if mod(max_steps,1)>0
  error('Please choose integer MAX_STEPS!')
end
  




%Extract relevent variables
[sensors,points] = size(data);
source = net.sources;
sources = length(source);
H_mean = net.posts.mix.H_mean;
H_prec = net.posts.mix.H_prec;
if POSITIVE
  [H,H_sq] = rect_expect(H_mean,H_prec);
else
  [H,H_sq] = deal(H_mean,H_mean.^2+1./H_prec);
end
lambda_hat = net.noise;
mn = net.mean;
dc = repmat(mn,1,points);
data - data-dc;
mn = net.mean;
mb = net.posts.noise.mb;
dc = repmat(mn,1,points);
data = data-dc;







%------------------------------Initialise------------------------------
% Initialise source signals
if init_x
  if POSITIVE
    old_x = rand(sources,points);
  else
    old_x = inv(H'*H)*H'*data;
  end
end

noiseH=repmat(lambda_hat,1,sources).*H;
data_bit = noiseH'*data;

preweight = (lambda_hat'*H_sq)';
weight_H = repmat(preweight,1,points);

for i = 1:sources
  betas{i} = source{i}.precs';
  mus{i} = source{i}.centres';
  comps(i) = source{i}.m;
  
  % Initialise source posterior precision
  beta_x_giv_q{i} = repmat(betas{i},1,points)+repmat(weight_H(i,:),comps(i),1);

  % Initialise source posterior means
  H_not_i = H;
  H_not_i(:,i) = 0;
  prior_mean = repmat(mus{i}.*betas{i},1,points);
  data_mean = data_bit(i,:) -  noiseH(:,i)'*(H_not_i*old_x);
  bracketbit = prior_mean + repmat(data_mean,comps(i),1);
  mu_x_giv_q{i} = 1./beta_x_giv_q{i}.*bracketbit;
  m = mu_x_giv_q{i};
  b = beta_x_giv_q{i};
  if POSITIVE
    [x_q{i},x_sq_q{i}] = rect_expect(m,b);
  else  
    [x_q{i},x_sq_q{i}] = deal(m,m.^2+1./b);
  end
  clear m b

  %work out hidden source states again
  source{i}.gammas = gammas(source{i},x_q{i},x_sq_q{i},2);
  x(i,:) = sum(source{i}.gammas.*x_q{i});
end  
%------------------------------Initialise------------------------------


oldold_x = x;
old_x = x;

err = tol+1;
steps = 0;
%--------------------------------unmix---------------------------------
while err>tol
  
  
  steps = steps+1;
  fprintf('%i/%i\n',steps,max_steps)
  
  %work out source posterior
  [x,x_sq,x_q,x_sq_q,mu_q,beta_q] = recon_source2(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps,POSITIVE,eta);

  %work out source state posterior probabilities
  for i=1:sources
    x_prec_q{i} = 1./(x_sq_q{i}-(x_q{i}).^2);
    source{i}.gammas = gammas(source{i},mu_q{i},beta_q{i},2);
    preentq(i,:) = -sum(source{i}.gammas.*log(source{i}.gammas+eps));
  end
  
 
  
  err = norm(old_x(:)-x(:));
  old_x = x;
  
  if steps == max_steps
    break;
  end



end
recon.signals = x;
recon.precs = 1./(x_sq-x.^2+eps);
recon.means_q = mu_q;
recon.precs_q = beta_q;
fprintf('State paths...')
for i=1:sources
  recon.state_probs{i} = source{i}.gammas;
  if HMM
    [block,ll,llb] = hmmdecode(x(i,:),x_sq(i,:),source{i});
    state_path(i,:) = block.q_star; %HMM state paths
  else
    [p,state_path(i,:)] = max(recon.state_probs{i});
  end
end
recon.state_paths = state_path;
fprintf('\n')
%--------------------------------unmix---------------------------------
