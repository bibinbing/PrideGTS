function log_like = ICA_like(x,hmm); 
% log_like = ICA_like(x,hmm);
%
% Finds average log likelihood of x under given vbHMM model.
%
% Called from 'hmmtrain'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% x           No. rows = No. sensors, No. columns = No. data points
% hmm         vbHMM-ICA model
%
%
%
% -----------
% Output
% -----------
%
% log_like    Average log likelihood of x under hmm
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
nets = hmm.obs.ica;
K = hmm.K;

for i=1:K
    pre_log_like(i,:) = avge_like(x,nets{i});
end

log_like = scalexp(pre_log_like)';

