function [b_alpha_hat,c_alpha_hat,f_alph] = update_alpha(priors,H_sq,src_type,ARD)
% [b_alpha_hat,c_alpha_hat,f_alph] = update_alpha(priors,H_sq,src_type,ARD)
%
% Works out ARD parameters for vbICA. 
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
% priors      Priors for ARD coefficients
% H_sq        Expecation of (mixmatrix.^2)
% src_type    Type of source model
% ARD:        1: Automatic Relevance Determination is ON  
%             0: Automatic Relevance Determination is OFF
%
%
%
% -----------
% Output
% -----------
%
% b,c _alpha_hat Posterior hyper-parameters
% f_alph         Contribution to free energy
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

[sensors sources nets] = size(H_sq);
if src_type=='e'
  pdf_fact = 1;
else
  pdf_fact=0.5;
end

b_alpha_0 = priors.b_alpha_0;
c_alpha_0 = priors.c_alpha_0;


%------------------------ARD coefficients alpha---------------------------


sum_mix_sq = sum(H_sq)';



if ARD
  %update hyper-parameters
  b_alpha_hat = 1./((1./b_alpha_0) + pdf_fact*sum_mix_sq);
  c_alpha_hat = c_alpha_0+sensors*pdf_fact;
  
  %calculate contribution to free energy.
  fa1 = sources*(gammaln(c_alpha_hat)-gammaln(c_alpha_0));
else
  %update hyper-parameters
  total_sq=sum(sum_mix_sq)*ones(sources,1);
  b_alpha_hat = 1./((1./b_alpha_0) + 0.5*total_sq);
  c_alpha_hat = c_alpha_0+sensors*sources/2;
  
  %calculate contribution to free energy.
  fa1 = sum(gammaln(c_alpha_hat)-gammaln(c_alpha_0));
end

fa2 = sum(c_alpha_hat*log(b_alpha_hat+eps)-c_alpha_0*log(b_alpha_0+eps));
f_alph = fa1+fa2;
%------------------------ARD coefficients alpha---------------------------
