function recon = unmix1(data,net,tol,max_steps,POSITIVE);
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

if tol<eps
  error('Please choose positive TOL!')
end

if max_steps<1
  error('Please choose MAX_STEPS > 1!')
end

if mod(max_steps,1)>0
  error('Please choose integer MAX_STEPS!')
end
%-------------------Error catching and default setting---------------------





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





%-----------------------------Initialise-----------------------------------
% Initialise source signals
if POSITIVE
  x = rand(sources,points);
else
  x = inv(H'*H)*H'*data;
end


for i = 1:sources
  betas{i} = source{i}.precs';
  comps(i) = source{i}.m;
  
  % Initialise hidden source states randomly
  gamms = rand(comps(i),points);
  gamms = gamms./repmat(sum(gamms),comps(i),1);
  
  % Initialise source posterior precision
  preweight = (lambda_hat'*H_sq)';
  weight_H = repmat(preweight,1,points);
  beta_x_hat(i,:) = (betas{i}'*gamms)+weight_H(i,:);
  
  if POSITIVE
    [x(i,:),x_sq(i,:)] = rect_expect(x,beta_x_hat);
  else  
    x_sq(i,:) = x(i,:).^2+1./beta_x_hat(i,:);
  end
  
  %work out hidden source states again
  source{i}.gammas = gammas(source{i},repmat(x(i,:),comps(i),1),repmat(x_sq(i,:),comps(i),1),1);
  
end  
%-----------------------------Initialise-----------------------------------



oldold_x = x;
old_x = x;

err = tol+1;
steps = 0;



%-------------------------------unmix--------------------------------------
while err>tol
  
  
  steps = steps+1;
  fprintf('%i/%i\n',steps,max_steps)
    
  %work out source posterior
  [x,x_sq,mu_x_hat,beta_x_hat,entx] = recon_source1(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps);

  %work out hidden source states
  for i=1:sources
    [source{i}.gammas,xi] = gammas(source{i},repmat(x(i,:),comps(i),1),repmat(x_sq(i,:),comps(i),1),1);
  end


  err = norm(old_x(:)-x(:));
  old_x = x;
  
  if steps == max_steps
    break;
  end



end


recon.signals = x;
recon.means = mu_x_hat;
recon.precs = beta_x_hat;
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
%-------------------------------unmix--------------------------------------






