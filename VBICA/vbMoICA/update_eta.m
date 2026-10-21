function [eta,ldl] = update_eta(data,nets,rho,x,x_prec,LEARN)
% [eta,ldl] = update_eta(data,nets,rho,x,x_prec,LEARN)
%
% This is the E-step of vbMoICA. 
%
% Called from 'vbmoica', 'unmix'. 
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% data      No. rows = No. sensors, No. columns = No. data points
% nets      Cell array containing the component ICA networks
% rho       Mixing prior hyperparameter after M-step
%
%
% Optional parameters
%
% x         Cell array of inferred source signals
% x_prec    Cell array of source signal precisions
% LEARN     1: Part of vbMoICA learning (Default)
%           0: Part of data unmixing
%
%
%
% -----------
% Output
% -----------
%
% eta       Calculated component responsibilities.
%
% Optional output
%
% ldl       Average log-likelihood per data point.

[sensors points] = size(data);
comps = length(rho);

if comps==1
  eta = ones(1,points);
  return
end


if nargin<5
    LEARN = 1;
else
    LEARN = 0;
end


for i=1:comps
  log_kappa_tild(i) = psi(rho(i))-psi(sum(rho));
  if LEARN
      like(i,:) = avge_like(data,nets{i});
  else
      like(i,:) = avge_like(data,nets{i},x{i},x_prec{i});
  end
end

big_lt_kappa = repmat(log_kappa_tild',1,points);
log_eta = (big_lt_kappa + like);


%scale to allow for machine precision.
eta = scalexp(log_eta);

if nargin>1
  ldl = like;
end

