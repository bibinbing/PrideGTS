function network = vbica1(data,states,source_type,SEE_PDFS,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX)
% network = vbica1(data,states,source_type,SEE_PDFS,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX)
% 
% ICA model with Automatic Relevance Determination (ARD) trained using 
% Variational Bayes with *Gaussian* source posterior.
% 
% y = Hx + e
%
% where each x has a MoG prior density and e is Gaussian noise.
% ARD scores each source on relevance and suppresses unnecessary sources.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        No. rows = No. sensors, No. columns = No. data points
% states      No. states per source model eg [2 4 3 5] for 4 sources
%
%
% Optional parameters: 
%
% source_type Mixture of 'g'aussians     (Default)
%             Mixture of 'e'xponentials
%             Mixture of 't'runcated gaussians
%             'hg','he','ht' - HMM versions
% SEE_PDFS    1: Observe evolution of source probability density functions
%             0: Don't observe evolution (Default)
% tol         Convergence tolerance.     (Default = 0.0001)   
% max_steps   max no_steps to cease.    (Default = 500)
% isonoise    1: learn isotropic noise.  (Default)
%             0: learn non-isotropic diagonal noise
% ARD:        1: Automatic Relevance Determination is ON   (Default)
%             0: Automatic Relevance Determination is OFF
% POS_SRC     1: learn positive bases. (Default for 't', 'e')
%             0: learn any bases.      (Default for 'g')
% POS_MIX     1: learn positive mixing matrix. (Default for 't', 'e')
%             0: learn any mixing matrix.      (Default for 'g')
%
%
%
% -----------
% Output
% -----------
%
% NETWORK is a data structure with the following fields:
%
% source_type   The source model
% source_states The source state vector
% sources       The source models
% recon         The reconstructed source signals (estimate of x)
% mixmatrix     The mixmatrix H
% alphas        The ARD parameters
% noise         The noise precision (precision of e)
% mean          The initialised noise mean
% priors        The model priors
% posts         The model posteriors
% state_paths   The hidden source state paths
% energy_path   The total NFE per iteration
% f_params      The NFE for the parameters
% f_hidd        The NFE for the hidden variables
% energy        The total negative free energy for the model
% model         'vbICA1'
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
  
%--------------------Error catching and default setting--------------------

addpath ../utils ../Shared ../MixModel1d
  
if nargin<2
  error('Please pass all necessary parameters!')
end

if nargin<3
  source_type = 'g';
end


global CHECK_PROGRESS;
if nargin<4
  CHECK_PROGRESS = 0;
else
  CHECK_PROGRESS= SEE_PDFS;
end

if nargin<5
  tol = 0.0001;
end

if nargin<6
  max_steps = 500;
end

if nargin<7
  isonoise = 1;
end

if nargin<8
  ARD = 1;
end

if length(source_type)==2
  src_type = source_type(2);
else
  src_type = source_type;
end

if nargin<9
  switch src_type
   case 'g'
    POS_SRC = 0;
   otherwise
    POS_SRC = 1;
  end
end

if nargin<10
  switch src_type
   case 'g'
    POS_MIX = 0;
   otherwise
    POS_MIX = 1;
  end
end


if sum(src_type == ['g' 't' 'e'])<1
  error('Please choose from available source_types!')
end


if tol<eps
  error('Please choose positive TOL!')
end

if max_steps<1
  error('Please choose MAX_STEPS > 1!')
end

if mod(max_steps,1)>0
  error('Please choose integer MAX_STEPS!')
end

if sum(isonoise == [0 1])~=1
  error('ISONOISE can only be ON (1) or OFF (0)!')
end

if sum(ARD == [0 1])~=1
  error('ARD can only be ON (1) or OFF (0)!')
end

if sum(POS_SRC == [0 1])~=1
  error('POS_SRC can only be ON (1) or OFF (0)!')
end

if sum(POS_MIX == [0 1])~=1
  error('POS_MIX can only be ON (1) or OFF (0)!')
end
%--------------------Error catching and default setting--------------------




%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
network = initialise_ica(data,source_type,states,isonoise);
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
network = learn_ica1(data,network,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX);
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::