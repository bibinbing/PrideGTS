function [gammas,xi] = gammas(src,y,y_sq,ALGORITHM)
% [gammas,xi] = gammas(src,y,y_sq,ALGORITHM)
%
% This is the E-step of mixmodel1d and the vbICA 
% source models. 
%
% Called from 'learn_mix1d', 'unmix1', 'unmix2'. 
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% src       Current model
% y         mixmodel1d: The observation data (sent as mxN)
%           vbICA1: Expectation of source signals (mxN) 
%           vbICA2: Location parameters of source posterior
% y_sq      mixmodel1d: y.^2 (sent as mxN)
%           vbICA: Expectation of (source signal.^2) (mxN)
%           vbICA2: Precisions of source posterior
% ALGORITHM 1: vbICA1  
%           2: vbICA2
%
%
%
% -----------
% Output
% -----------
%
% gammas    Calculated posterior probabilties
% xi        Calculated trans. probs for HMM  (otherwise = 0)
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

  
littlebit = eps;
HMM=0;

  
%------------------------Extract appropriate variables---------------------
[m,N] = size(y);
src_type = src.type;
lambda = src.posts.lambda;
lambda_tot=sum(lambda);
%------------------------Extract appropriate variables---------------------



%----------------------------Data pseudo-likelihood------------------------
if ALGORITHM==1  
  % mixmodel1d and vbICA1
  log_obslike_tilde = log_ptilde1(y,y_sq,src);
else             
  % vbICA2
  log_obslike_tilde = log_ptilde2(y,y_sq,src);
end
%----------------------------Data pseudo-likelihood------------------------



%------------------------------Extra bit for HMM---------------------------
if length(src_type) == 2
  HMM=1;
  eta = src.posts.eta;
  etasum = sum(eta,2);
  % 'Tilded' transition matrix for FB part
  for i=1:m
      for j=1:m
	  index = digamma(eta(i,j))-digamma(etasum(i));
	  P_tilde(i,j) = exp(index);
      end
  end
end
%------------------------------Extra bit for HMM---------------------------



%----------------------------Work out probabilities------------------------
for s=1:m,
   log_tilde_pi(s)=psi(lambda(s))-psi(lambda_tot);
end

if HMM
  % Forward-Backward bit if HMM
  [gammas,xi] = FB(exp(log_tilde_pi),(scalexp(log_obslike_tilde))',P_tilde);
else
  log_pi_tilde = repmat(log_tilde_pi',1,N);
  log_gam = log_pi_tilde + log_obslike_tilde;
  if m==1
    gammas = scalexp(log_gam,1);  % scale to get past machine precision
  else
    gammas = scalexp(log_gam);
  end
  xi = 0;
end
%----------------------------Work out probabilities------------------------
