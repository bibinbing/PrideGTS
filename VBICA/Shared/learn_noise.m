function [b_Lam_hat,c_Lam_hat,f_lambda] = learn_noise(priors,data,mb,H,x,H_sq,x_sq,isonoise,eta)
% [lambda_hat,b_Lam_hat,c_Lam_hat,f_lambda] = learn_noise(priors,data,H,x,H_sq,x_sq,isonoise,eta)
%
% Works out noise precision posterior for vbICA. 
%
% Called from 'learn_ica1', 'learn_ica2.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% priors      Priors for noise precision
% data        No. rows = No. sensors, No. columns = No. data points
% mb          Precision of noise mean posterior
% H           Expecation of mixmatrix
% x           Expectation of source signals
% H_sq        Expecation of (mixmatrix.^2)
% x_sq        Expectation of (source signals.^2)
%
%
%
% Optional parameters
%
% isonoise    1: learn isotropic noise.  (Default)
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
% b,c _Lam_hat  Parameters for noise precision posterior
% f_lambda      Contribution to free energy
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
if nargin<8
  isonoise = 1;
end

if nargin<9
  eta = 1;
end


[sensors points] = size(data);
[sensors sources] = size(H);

b_Lam_0 = priors.b_Lam_0;
c_Lam_0 = priors.c_Lam_0;

sumeta = points;

if length(eta)~=1
  sumeta = sum(eta);
  eta_sens = repmat(eta,sensors,1);
else
  eta_sens = 1;
  sumeta = points;
end



%----------------------Sensor noise precision Lambda-----------------------
data_hat = H*x;
dd=eta_sens.*data.^2;
ddh=eta_sens.*(data.*data_hat);

data_hat_var = H_sq*x_sq-(H.^2)*(x.^2);
dhdh = eta_sens.*(data_hat.^2+data_hat_var);

mean_bit = repmat(1./mb,1,points).*eta_sens;
Ed = sum(dd-2*ddh+dhdh+mean_bit,2);


if (isonoise)
 
  Ed = sum(Ed);

  %update hyper-parameters
  b_Lam_hat =1./(1/b_Lam_0 + 0.5*Ed);
  c_Lam_hat = c_Lam_0+(sumeta*sensors/2);
  
  %calculate contribution to free energy.
  fl1=gammaln(c_Lam_hat)-gammaln(c_Lam_0);
  fl2=c_Lam_hat*log(b_Lam_hat)-c_Lam_0*log(b_Lam_0);
  f_lambda=fl1+fl2;

  b_Lam_hat = repmat(b_Lam_hat,sensors,1);
else
  %update hyper-parameters
  b_Lam_hat =1./(1/b_Lam_0 + 0.5*Ed);
  c_Lam_hat = c_Lam_0+(sumeta/2);
  
  %calculate contribution to free energy.
  fl1=gammaln(c_Lam_hat)-gammaln(c_Lam_0);
  fl2=c_Lam_hat.*log(b_Lam_hat)-c_Lam_0.*log(b_Lam_0);
  f_lambda=sum(fl1+fl2);
end

clear fl1 fl2;
%----------------------Sensor noise precision Lambda-----------------------