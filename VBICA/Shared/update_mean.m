function [mn,mb,f_dc] = update_mean(priors,data,H,x,lambda_hat,eta)
% [mn,mb,f_dc] = update_mean(priors,data,H,x,lambda_hat,eta)
%
% Works out noise mean posterior for vbICA. 
%
% Called from 'learn_ica1', 'learn_ica2'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% priors      Priors for noise mean
% data        No. rows = No. sensors, No. columns = No. data points
% H           Expecation of mixmatrix
% x           Expectation of source signals
% lambda_hat  Expectations of noise precision
%
%
%
% Optional parameters
%
% eta         Class responsibilities if being trained as part of a 
%             Mixture of ICAs, otherwise = 1
%
%
%
% -----------
% Input
% -----------
%
% mn,mb       Posterior parameters
% f_dc        Contribution to free energy
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)  

if nargin<6
  eta = 1;
end

[sensors points] = size(data);

if length(eta) == 1
  sumeta = points;
else
  sumeta = sum(eta);
  eta = repmat(eta,sensors,1);
end

mb0 = priors.mb0;
mn0 = priors.mn0;



%work out precision of noise mean posterior
mb = mb0 + sumeta*lambda_hat;

%work out mean of noise mean posterior
ratio = lambda_hat./mb;
diff = mb0.*mn0+sum(eta.*(data - H*x),2);
mn = ratio.*diff;

%contribution to free energy
pratio = mb0./mb;
bit1 = pratio-1-log(pratio);
bit2 = mb0.*(mn-mn0).^2;
f_dc = -0.5*sum(bit1+bit2);