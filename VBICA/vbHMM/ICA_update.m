function hmm = ICA_update(x,gamma,hmm)
% hmm = ICA_update(x,gamma,hmm);
%
% Learn the constituent ICA networks of hmm.
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
% gamma       Posterior state probabilities
% hmm         vbHMM-ICA model
%
%
%
% -----------
% Output
% -----------
%
% hmm         Updated HMM data structure
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)    
  
global CHECK_PROGRESS;

nets = hmm.obs.ica;
tol = hmm.train.tol;
max_steps = hmm.train.cyc;
isonoise = hmm.obs.isonoise;
ARD = hmm.obs.ARD;
POS_SRC = hmm.obs.POS_SRC;
POS_MIX = hmm.obs.POS_MIX;

[Fica nets] = learn_mica(x,nets,gamma,tol,max_steps,ARD,isonoise,POS_SRC,POS_MIX);

hmm.obs.ica = nets;
hmm.obs.energy = sum(Fica);