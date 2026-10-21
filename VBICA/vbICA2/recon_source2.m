function [x,x_sq,x_q,x_sq_q,mu_x_giv_q,beta_x_giv_q,entx] = recon_source2(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps,POSITIVE,eta);
% [x,x_sq,x_q,x_sq_q,mu_x_giv_q,beta_x_giv_q,entx] = recon_source2(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps,POSITIVE,eta);
%
% Works out source posterior (Mixture of Gaussians) for vbICA2. 
%  
% Called from 'learn_ica2', 'unmix2' of vbICA2.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        No. rows = No. sensors, No. columns = No. data points
% old_x       Previous expecation of hidden source posterior
%             Used as starting point for new calculation
%             No. rows = No. sources, No. columns = No. points
% H,H_sq      Expectations of mix matrix
% lambda_hat  Expectation of noise precision
% source      Cell array of source models
%
%
% Optional parameters:
%
% tol         Convergence tolerance (Default = 1e-4)
% max_steps   Max no_steps to cease (Default = 500) 
% POSITIVE    Positive-only source signals (Default = 0 for 'g')
%                                          (Default = 1 for 'e', 't')
% eta         Class responsibilities if being trained as part of a 
%             Mixture of ICAs. Otherwise = 1;
%
%
% -----------
% Output
% -----------
%
% x,x_sq       New expectations of hidden source vector under posterior
% x_q,x_sq_q   New expectations of hidden source vector under posterior 
%              given source states
% mu_x_giv_q   Location parameters of posterior
% beta_x_giv_q Precisions of source posterior
% entx         Contribution to free energy
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

%--------------------Error catching and default setting--------------------
if nargin<6
  error('Please pass all necessary parameters!')
end

if nargin<7
  tol = 0.0001;
end

if nargin<8
  max_steps = 500;
end

source_type = source{1}.type;
if length(source_type)==2
  src_type = source_type(2);
else
  src_type = source_type;
end
if nargin<9
  switch src_type
   case 'g'
    POSITIVE = 0;
   otherwise
    POSITIVE = 1;
  end
end

if nargin<10
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
%--------------------Error catching and default setting--------------------  
[sensors sources] = size(H);
[sensors points] = size(data);
MIXTURE = 0;

if length(eta) ~= 1
  eta = repmat(eta,sources,1);
  MIXTURE = 1;
end


%--------------------------------P(x)--------------------------------------


%Work out beta_hat_i,q (p'(x|q))
preweight = (lambda_hat'*H_sq)';
weight_H = repmat(preweight,1,points).*eta;
for i = 1:sources
  comps(i) = source{i}.m;
  gammas{i} = source{i}.gammas;

  switch src_type
   case 'g'
    betas{i} = source{i}.precs';
    mus{i} = source{i}.centres';
    mean_prior{i} = betas{i}.*mus{i};
   case 'e'
    mean_prior{i} = -source{i}.precs';
    betas{i} = zeros(comps(i),1);
   case 't' 
    betas{i} = source{i}.precs';
    mean_prior{i} = zeros(comps(i),1);
  end

  %work out precision of P(x_i|q_(i,q_i),data)
  beta_x_giv_q{i} = repmat(betas{i},1,points)+repmat(weight_H(i,:),comps(i),1);
end
clear preweight weight_H 
  
 
%work out mu_hat_i,q (p'(x_q))
noiseH=repmat(lambda_hat,1,sources).*H;
data_bit = noiseH'*data;

oldold_x = old_x(:);
errX = tol+1;
xsteps = 0;

while errX>tol
  %fprintf('x');
  xsteps = xsteps+1;
  clear new_x x;
   
  for i=1:sources

    if ~MIXTURE
      eta(i,1:points) = 1;
    end
    
    
    H_not_i = H;
    H_not_i(:,i) = 0;
    prior_mean = repmat(mean_prior{i},1,points);
    data_mean = data_bit(i,:) -  noiseH(:,i)'*(H_not_i*old_x);
    bracketbit = prior_mean + repmat(eta(i,:).*data_mean,comps(i),1);
    
    
    %work out mean of P(x_i|q_(i,q_i),data)
    mu_x_giv_q{i} = 1./beta_x_giv_q{i}.*bracketbit;
    
    clear prior_mean bracketbit
    
    if POSITIVE
      x_q{i} = rect_expect(mu_x_giv_q{i},beta_x_giv_q{i});
    else
      x_q{i} = mu_x_giv_q{i};
    end
    
    clear m b
    
    
    %work out mean of P(x_i|data)
    x(i,:) = sum(gammas{i}.*x_q{i},1);
   
    
  end
 
  errX=norm(x(:)-old_x(:));
  
  if xsteps==max_steps
    break;
  end
  old_x = x;
end




%free energy contribution
if POSITIVE
  for i=1:sources
    [x_q{i} x_sq_q{i}] = rect_expect(mu_x_giv_q{i},beta_x_giv_q{i});
    x_sq(i,:) = sum(gammas{i}.*x_sq_q{i},1);
    prepreentx = NFEpos(mu_x_giv_q{i},beta_x_giv_q{i});
    preentx(i,:) = sum(gammas{i}.*prepreentx);
  end
  entx = sum(sum(preentx));
else
  for i=1:sources
    x_sq_q{i} = mu_x_giv_q{i}.^2+1./beta_x_giv_q{i};
    x_sq(i,:) = sum(gammas{i}.*x_sq_q{i},1);
    preentx(i,:) = -sum(gammas{i}.*log(beta_x_giv_q{i}),1);
  end
  entx = 0.5*sum(sum(preentx,1))+sources*points/2;
end

 
%--------------------------------P(x)--------------------------------------
