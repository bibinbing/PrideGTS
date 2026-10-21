function probs = plotMoG(m,s,p,xmin,xmax,normfact);
% probs = plotMoG(m,s,p,xmin,xmax,normfact);
%
% Plots a 1d Mixture of Gaussians.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% m        MIX1 model or Vector of means
% s        Vector of variances
% p        Vector of mix coefficients
%
% NOTE     If a Mix1 data structure is passed, then s is
%          is an optional parameter representing minimum
%          x-axis value and p an optional max value
%
%
% Optional parameters
%
% xmin     Minimum x-axis value
% xmax     Maximum x-axis value
% Normfact Adjustment when used for rectified Gaussians
%
%
%
% -----------
% Output
% -----------
%
% Optional output
%
% probs Probabilities calculated for the overall curve
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
SETWIDTH = 1;
    
if isstruct(m)
  mog = m;
  m = mog.centres;
  if nargin>1
      xmin = s;
      xmax = p;
      SETWIDTH = 0;
  else
      SETWIDTH = 1;
  end
  s = 1./mog.precs;
  p = mog.pi;
end


if nargin<6
  normfact = 1;
end

if sum(p)~=1
  p = p./sum(p);
end

if SETWIDTH
    if nargin<4
	width=7;
	dm = max(m) - min(m) + width*p*sqrt(s');
	if dm<max(s)
	    dm = 5*p*sqrt(s');
  end
  mm = mean(m);
  xmin = mm - dm;
  xmax = mm + dm;
    end
end

dt = (xmax-xmin)/1000;

x=[xmin:dt:xmax];


comps = length(p);
prob=0;

for n=1:comps
  norming = normfact/(2*pi*s(n)).^0.5;
  expfactor = exp(-0.5*(1/s(n))*(x-m(n)).^2);
  compprob(:,n+1) = (p(n)*norming.*expfactor)';
  prob = prob+compprob(:,n+1);
end


compprob(:,1)=prob;
plot(x,compprob(:,2:comps+1),'--','linewidth',2);
hold on
plot(x,compprob(:,1),'k-','linewidth',2);
hold off

if nargout>0
  probs = compprob(:,1);
end
