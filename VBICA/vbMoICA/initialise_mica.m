function [eta,nets] = initialise_mica(data,states,source_type,isonoise);
% [eta,nets] = initialise_mica(data,source_type,states,isonoise);
% 
% Intialisation for ICA mixture model. A Gaussian mixture model
% is used to partition the data into length(states) clusters. 
% The associated probabilities are hard-limited to 1 or 0 and 
% fed to the intialise_ICA function.
%
% Called from 'vbmoica', 'hmminit' (vbHMM).
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        no_rows=no_sensors. no_columns=no_data points.
% states      Cell array that contains number states per ICA source MoG 
%             eg STATES{1} = [2 4 3 5] means ICA model 1 has 4 sources with
%             2,4,3,5 components respectively
% source_type Mixture of 'g'aussians     (Default)
%             Mixture of 'e'xponentials
%             Mixture of 't'runcated gaussians
%             'hg','he','ht' - HMM versions
%
% Optional prameters:
%
% isonoise    1: learn isotropic noise.  (Default)
%             0: learn non-isotropic diagonal noise
%
%
%
% -----------
% Output 
% -----------
%
% eta         The hard-limited partition probabilities
% nets        The cell array of individual ICA models   
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
global CHECK_PROGRESS; 

if nargin < 4
  isonoise = 1;
end
fprintf('Initialising mixture probabilities using a GMM...\n\n')
comps = length(states);
[sensors points] = size(data);
%src_type = source_type(length(source_type));
%if src_type == 'g'
  warning = 1;
  % Initilaise mix probs using Gaussian mix model.
  % Needs to be improved for vbMoICA with positive sources.
  while warning
    warning = 0;
    [mix,m,c,p,ev] = my_gmm(data',comps,'f');
    warning = mix.warning;
    if warning
      fprintf('GMM initialisation did not converge - Retrying...')
    end
  end

%else
%  eta = rand(comps,points);
%  sumeta = repmat(sum(eta,1),comps,1);
%  eta = eta./sumeta;
%end

  eta = mix.probs'+eps;
 
  init_eta = round(eta==repmat(max(eta),comps,1));
  
  % Initialise individual ICA networks
  for i=1:comps
      if CHECK_PROGRESS
	  figure(gcf+1);
      end
      nets{i} = initialise_ica(data,source_type,states{i},isonoise,init_eta(i,:));
  end
