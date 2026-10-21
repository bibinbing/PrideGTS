function [hmm,Fp,Fpi,entHidden]=hmm_mstep(hmm,Gamma,Xi);
% [hmm,Fp,Fpi,entHidden]=hmm_mstep(hmm,Gamma,Xi);
%
% Learns the HMM-specific parameters in a 1-dimensional HMM 
% using the variational Bayes framework.
%
% Called from 'learn_mix1d'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% hmm        Current model
% Gamma      State probabilities
% Xi         Pairwise transition probabilities
%
%
%
% -----------
% Output
% -----------
%
% hmm        HMM model. P: transition matrix
%                       pi: initial state probabilities
%                       posts.eta: posterior params for P
%                       posts.lambda: posterior params for pi
% Fp         Contribution to energy from trans. matrix
% Fpi        Contribution to energy from initial state probs.
% entHidden  Entropy of hidden states
%
%
% --------------------------------------------------------------
%
% Original code by Iead Rezek
% Paper: Learning Ensemble Hidden Markov Models for Biosignal 
%        Analysis (www.robots.ox.ac.uk/~parg)
%
% Modified by Rizwan Choudrey for use in vbICA model
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
  [K,T] = size(Gamma);  
  
  %-----------------------Transition matrix------------------------------
  sxi=sum(Xi);
  sxi=reshape(sxi,K,K);
  eta0 = hmm.priors.eta_0;
  eta = eta0+sxi;                                  %'iota_{cd}' in thesis
  sumeta = repmat(sum(eta,2),1,K);
  trueP=eta./sumeta;                               %'R' in thesis
    
  %contribution to free energy
  if length(eta0(1,:))~=K
      eta0 = repmat(eta0,1,K);
  end
  bit1 = sum(gammaln(eta) - gammaln(eta0),2);
  bit2 = gammaln(sum(eta,2))-gammaln(sum(eta0,2));
  Fp = sum(bit1-bit2);
  %-----------------------Transition matrix------------------------------
  
  
  
  %----------------------Initial state probs-----------------------------
  lambda0 = hmm.priors.lambda_0;
  lambda = lambda0+Gamma(:,1)';
  sumlambda = sum(lambda);
  truePi = lambda./sum(lambda);
     
  %contribution to free energy
  if length(lambda0)~=K
      lambda0 = repmat(lambda0,1,K);
  end
  bit1 = sum(gammaln(lambda) - gammaln(lambda0));
  bit2 = gammaln(sum(lambda))-gammaln(sum(lambda0));
  Fpi = sum(bit1-bit2);
  %----------------------Initial state probs-----------------------------
  
  
  
  %---------------------Contribution to Energy---------------------------
  entStart = -sum(Gamma(:,1).*log(Gamma(:,1)+eps));
  for i=1:K
    Xi3d = reshape(Xi,[T-1,K,K]);
    sxi = sum(Xi3d(:,:,i),2);
    Psi(:,:,i) = Xi3d(:,:,i)./(repmat(sxi,1,K)+eps);
  end
  entXi = -sum(sum(sum(Xi3d.*log(Psi+eps))));
  entHidden = entStart+entXi;
  %---------------------Contribution to Energy---------------------------
  
  
  
  hmm.P = trueP;
  hmm.pi = truePi;
  hmm.posts.eta = eta;
  hmm.posts.lambda = lambda;