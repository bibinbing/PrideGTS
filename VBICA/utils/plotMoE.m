function probs = plotMoE(b,p,xmax);
% probs = plotMoE(b,p,xmax) 
%
% Plots a 1d Mixture of Exponentials.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% b     MIX1 model or vector of scale parameters
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
if isstruct(b)
  mog = b;
  if nargin>1
    xmin = p;
    SETWIDTH = 0;
  else
    SETWIDTH = 1;
  end
  b = 1./mog.precs;
  p = mog.pi;
end

p = p./sum(p);

  
if nargin<3
  width=2;
  %xmax=200;
  xmax = width*max(b).^2;
end


dt = (xmax-xmin)/1000;

x=[xmin:dt:xmax];


comps = length(p);
prob=0;

for n=1:comps
  norming = 1/b(n);
  expfactor = exp(-x./b(n));
  compprob(:,n+1) = (p(n)*norming.*expfactor)';
  prob = prob+compprob(:,n+1);
end

compprob(:,1)=prob;
plot(x,compprob(:,2:comps+1),'--');
hold ("on")
plot(x,compprob(:,1),'k-','linewidt',2);
hold ("off")
if nargout>0
  probs = compprob(:,1);
end