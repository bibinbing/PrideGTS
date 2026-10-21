function [recons,eta] = unmix(data,net,tol,max_steps,POSITIVE);

% [recons eta] = unmix(data,net,tol,max_steps);
% 
% Use model trained by Mix_vbICA to unmix data and infer hidden
% source signals and states.
%
%
% -----------
% Input 
% -----------
%
% Necessary parameters
%
% data        No. rows = No. sensors, No. columns = No. data points
% net         MoICA network to use for unmixing
%
%
%
% Optional Parameters
%
% tol         Convergence tolerance  (Default=1e-4)  
% max_steps   Max no_steps to cease  (Default=100)
% POSITIVE    Positive-only source signals (default=1 for 'e', 't')
%                                          (default=0 for 'g')
%
%
%
% -----------
% Output 
% -----------
%
% recons      Cell array of RECON data structures
% eta         Posterior probabilities for Data under each ICA
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)  


if nargin<3
  tol = 1e-4;
end

if nargin<4
  max_steps = 100;
end

if nargin<5
  src_type = net.nets{1}.source_type(length(net.nets{1}.source_type));
  switch src_type
   case 'g'
    POSITIVE = 0;
   case 'e'
    POSITIVE = 1;
   case 't'
    POSITIVE = 1;
  end
end

[sensors points] = size(data);
nets = net.nets;
comps = length(nets);
rho = net.posts.rho;
rho_0 = net.priors.rho_0;
err = tol+1;


%-----------------------------Initialise-----------------------------------
steps = 1;
eta = ones(comps,points)./comps;
for i=1:comps
  fprintf('Reconstructing using ICA %i...\n',i)
  recons{i} = unmix2(data,nets{i},tol,max_steps,POSITIVE,[],eta(i,:));
  old_x{i} = recons{i}.signals;
  x_prec{i} = recons{i}.precs;
end

oldeta = eta(:);
fprintf('Updating mixture coefficients...\n\n')
[eta avge_like] = update_eta(data,nets,rho,old_x,x_prec);
err = norm(oldeta-eta(:))/length(eta);
%-----------------------------Initialise-----------------------------------



%-----------------------------Reconstruct----------------------------------
while err>tol
  steps = steps+1;
  %-----------------------------Mix nets--------------------------------
  for i=1:comps
    fprintf('Reconstructing using ICA %i...\n',i)
    recons{i} = unmix2(data,nets{i},tol,max_steps,POSITIVE,old_x{i},eta(i,:));
    x{i} = recons{i}.signals;
    x_prec{i} = recons{i}.precs;
    old_x{i} = x{i};
  end
  %-----------------------------Mix nets--------------------------------
  
  
  %-----------------------------Mix probs-------------------------------
  fprintf('Updating mixture coefficients...\n\n')
  oldeta = eta(:);
  [eta avge_like] = update_eta(data,nets,rho,x,x_prec,0);
  err = norm(oldeta-eta(:))/length(eta);
  fprintf('%i/%i.....Change in recon = %1.4f percent\n\n',steps,max_steps,err*100);
  %-----------------------------Mix probs-------------------------------

  
  if steps == max_steps
    break;
  end  
  
end
%----------------------------Reconstruct-----------------------------------

