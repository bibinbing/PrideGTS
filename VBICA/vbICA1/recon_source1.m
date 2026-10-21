function [x,x_sq,mu_x_hat,beta_x_hat,entx] = recon_source1(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps,POSITIVE);
% [x,x_sq,mu_x_hat,beta_x_hat,entx] =  recon_source1(data,old_x,H,H_sq,lambda_hat,source,tol,max_steps,POSITIVE); 
%
% Works out source posterior (Gaussian) for vbICA1. 
%  
% Called from 'learn_ica1', 'unmix1' of vbICA1.
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
% tol         Convergence tolerance  (Default = 1e-4)
% max_steps   Max no_steps to cease  (Default = 500) 
% POSITIVE    Positive-only source signals (Default = 0 for 'g')
%                                          (Default = 1 for 'e', 't')
%
%
%
% -----------
% Output
% -----------
%
% x,x_sq      New expectations of hidden source vector under posterior
% mu_x_hat    Location parameter of posterior
% beta_x_hat  Precision of posterior
% entx        Contribution to negative free energy
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

%--------------------------------p'(x)-------------------------------------
for i=1:length(source)
    gamma = source{i}.gammas;
    precision = source{i}.precs;    
    switch src_type
     case 'g'
      srcmean = source{i}.centres;
      beta_bar(i,:) = precision*gamma;
      mu_bar(i,:) = (precision.*srcmean)*gamma;
     case 'e'
      beta_bar = 0;
      mu_bar(i,:) = -precision*gamma;
     case 't' 
      beta_bar(i,:) = precision*gamma;
      mu_bar = 0;
    end
    
    clear gamma precision srcmean;
    
end

%Work out beta_hat_i (p'(x))
postsum = (lambda_hat'*H_sq)';
beta_x_hat = beta_bar + repmat(postsum,1,points);
clear postsum H_sq;


%work out mu_hat_i (p'(x))
noiseH=repmat(lambda_hat,1,sources).*H;
h1=diag(noiseH'*H);
for  iter=1:max_steps 
  clear new_x x;
  %recognising that y-y_hat(k!=i) = y-[y_hat(all i)-h_ji*x_i]
  h2 = noiseH'*(data-H*old_x);
  h3 = repmat(h1,1,points).*old_x + h2;
  bracketbit = mu_bar + h3;
  

  mu_x_hat = (1./beta_x_hat).*bracketbit;

  if POSITIVE
    x = rect_expect(mu_x_hat,beta_x_hat);
  else
    x = mu_x_hat;
  end

  old_x = old_x(:);
  new_x = x(:);
  errX=norm(new_x-old_x);

  if errX < tol
      break;
  end
  old_x = x;
end

if POSITIVE
  [x,x_sq] = rect_expect(mu_x_hat,beta_x_hat);
  preentx = NFEpos(mu_x_hat,beta_x_hat);
  entx = sum(sum(preentx));
else
  x_sq = mu_x_hat.^2+1./beta_x_hat;
  preentx = -0.5*log(beta_x_hat);
  entx = sum(sum(preentx,2))+points*sources/2; 
end


%---------------------------------p'(x)------------------------------------
