function probs = plotSrc(mix1);
% probs = plotSrc(mix1)
%
% Plots a 1d Mixture model learnt via MixModel1d.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% mix1   MixModel1d mixture model 
%
%
%
% -----------
% Output
% -----------
%
% Optional output
%
% probs  Probabilities calculated for the overall curve
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
if nargin < 1
  error('Please pass all necessary parameters!')
end

if ~isStruct(mix1)
  error('Please pass a MixModel1d data structure!')
end



source_type = mix1.type;
if length(source_type)==2
  HMM = 1;
  src_type = source_type(2);
  disp('HMM - plotting overall density using sum(gammas)...')
  lambda = sum(mix1.gammas');
  mix1.pi = lambda./sum(lambda);
else
  HMM = 0;
  src_type = source_type;
end

if nargout>1
  switch src_type
   case 'g'
    probs = plotMoG(mix1);
   case 'e'
    probs = plotMoE(mix1);
   case 't'
    probs = plotRMoG(mix1);
  end
else
  switch src_type
    case 'g'
    plotMoG(mix1);
   case 'e'
    plotMoE(mix1);
    case 't'
    plotRMoG(mix1);
   end
end

