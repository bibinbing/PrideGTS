function [mix,m,C,p,ev] = my_gmm(x,K,f_or_d);
% [mix,m,C,p,ev] = my_gmm(x,K,full_or_diag);
%
% Modification of NETLAB's Gaussian mixture model front-end.
% For use in intialising vbMoICA.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% x         M x T data matrix
% K         Number of Gaussian components
% f_or_d    'f'ull or 'd'iagonal covariances
%
%
% -----------
% Output
% -----------
% 
% MIX       NETLAB's GMM data structure
% m         Component means
% C         Component covariances
% p         Component mixture probabilities
% ev        Posterior component probabilities
%
%
% --------------------------------------------------------------
%
% Original code Copyright (c) Christopher M Bishop, Ian T Nabney
% (1996, 1997). Available in NETLAB distribution.
% www.ncrg.aston.ac.uk/netlab
%
% Modified by Rizwan Choudrey for use in vbICA model
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg) 
  
  
[N,D] = size(x);
t = zeros(N,1);

if ( (f_or_d == 'f')  | (f_or_d == 'd'))
else
  error('full_or_diag must be f or d');
  return;
end;

if (f_or_d == 'f')
  mix = gmm(D,K,'full');
else
  mix = gmm(D,K,'diag');
end;
   
options = foptions;
options(14) = 5;		% iterations of K-means
mix = gmminit(mix, x, options);
options(1)  = 0;		% Prints out error values.
options(5) = 0;		% trap bad covariances
options(14) = 100;		% Max. number of iterations.
mix = my_gmmem(mix, x, options);
pp = gmmpost(mix,x);
mix.probs=pp;

ev = gmmprob(mix,x);
m = mix.centres;
if (f_or_d == 'f')
  C=[];for n=1:K,C=[C;mix.covars(:,:,n)];end;
else
  C = mix.covars;
end;
p = mix.priors;
return;