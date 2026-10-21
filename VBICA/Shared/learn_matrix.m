function [H,H_sq,m_H_hat,alpha_hat,f_mix] = learn_matrix(data,old_H,x,x_sq,lambda_hat,alpha_prec,alpha_mean,tol,max_steps,POSITIVE,eta)
% [H,H_sq,m_H_hat,alpha_hat,f_mix] = learn_matrix(data,old_H,x,x_sq,lambda_hat,alpha_prec,alpha_mean,tol,max_steps,POSITIVE,eta)
%
% Works out mixmatrix posterior. Called from 'learn_ica1', 'learn_ica2.
%
%
% Necessary parameters
%
% data        No. rows = No. sensors, No. columns = No. data points
% old_H       Previous expecation of mixmatrix
%             Used as starting point for new calculation
% x,x_sq      Expectations of hidden source signals
% lambda_hat  Expectations of noise precision
% alpha_prec  Expectation of ARD coefficient (zero for 'e')
% alpha_mean  Expectation of ARD coefficient for 'e' (zero otherwise)
%
%
% Optional parameters
%
% tol         Convergence tolerance     (Default = 0.0001)   
% max_steps   max no_steps to cease     (Default = 500)
% POSITIVE    Positive-only source signals (Default = 0 for 'g')
%                                          (Default = 1 for 'e', 't')
% eta         Class responsibilities if being trained as part of a 
%             Mixture of ICAs, otherwise = 1
%
%
%
% -----------
% Output
% -----------
%
% H,H_sq      New expectations of mixmatrix under posterior
% m_H_hat     Location parameter of posterior
% alpha_hat   Precision of mixmatrix posterior
% fmix        Contribution to free energy
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

if nargin<8
  tol = 1e-4;
end

if nargin<9
  max_steps = 500;
end

if nargin<10
  POSITIVE = 0;
end

if nargin<11
  eta = 1;
end



[sensors points] = size(data);
[sensors sources] = size(old_H);

if length(eta)~=1
  eta = repmat(eta,sources,1);
end


%-----------------------------------P(H)------------------------------------
%work out alpha_hat, precision of posterior.
sum_x_sq = sum(eta.*x_sq,2)';
network_term = repmat(lambda_hat,1,sources).*repmat(sum_x_sq,sensors,1);
alpha_hat = repmat(alpha_prec',sensors,1) + network_term;
clear x_sq sum_x_sq network term;

%work out m_H_hat (p'(h))
coefficient = 1./alpha_hat;
big_noise = repmat(lambda_hat,1,sources);
x1 = diag(eta.*x*x');

for  iter=1:max_steps
  clear new_H H;
  
  %remembering that (y-y_hat(k!=i)) = [y-(y_hat(all_i)-h_ij*x_i)]
  x2 = (data-old_H*x)*(x.*eta)';
  x3 = repmat(x1',sensors,1).*old_H + x2;
  
  old_H = old_H(:);
  
  %work out mean of posterior.
  m_H_hat = coefficient.*(-repmat(alpha_mean',sensors,1)+big_noise.*x3);
  %  m_H_hat = coefficient.*(big_noise.*x3);
  
  if POSITIVE
    H = rect_expect(m_H_hat,alpha_hat);
  else
    H = m_H_hat;
  end
  
  new_H = H(:);

  errH=norm(new_H-old_H);
    
  if errH < tol
    break;
  end
  old_H = H;
end 

if POSITIVE
  [H,H_sq] = rect_expect(m_H_hat,alpha_hat);
  f_mix = sum(sum(NFEpos(m_H_hat,alpha_hat)));
else
  H_sq = m_H_hat.^2+1./alpha_hat ;
  fa4=-0.5*sum(sum(log(alpha_hat)));
  f_mix=0.5*sensors*sources+fa4;
end
%-----------------------------------P(H)----------------------------------------