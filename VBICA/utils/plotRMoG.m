function plotRMoG(s,p,xmax);
% probs = plotRMoG(s,p,xmax);
%
% Plots a 1d Mixture of rectified Gaussians.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% s     MIX1 model or vector of variance parameters
% p     Vector of mixing parameters
%
% NOTE  If a Mix1 data structure is passed, then p is
%       is an optional parameter representing minimum
%       x-axis value
%
%
% Optional parameters
%
% xmax  Maximum x-axis value
%
%
%
% -----------
% Output
% -----------
%
% Optional output:
%
% probs Probabilities calculated for the overall curve
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

xmin =0;
if isstruct(s)
  mog = s;
  if nargin>1
    xmin = p;
    SETWIDTH = 0;
  else
    SETWIDTH = 1;
  end
  s = 1./mog.precs;
  p = mog.pi;
end
  
m = zeros(1,length(s));

if nargin<3
  width=5;
  [maxm,k] = max(m);
  [maxs,j] = max(s);
  xmax = ceil(max([maxm+width*s(k),m(j)+width*maxs]));
end

plotMoG(m,s,p,xmin,xmax,2);