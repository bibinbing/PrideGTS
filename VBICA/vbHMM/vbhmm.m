function hmm = vbhmm(data,states,source_type,SEE_PDFS,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets);
% hmm = vbhmm(data,states,source_type,SEE_PDFS,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets)
%
% ICA-based Hidden Markov Model trained by Variational Bayes inc. 
% Automatic Relevance Determination (ARD).
% 
% ICA model: y = Hx + e
%
% where each x has a mixture density and e is isotropic
% Gaussian noise with precision lambda.
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
% states      Cell array that contains number states per ICA source MoG 
%             eg STATES{1} = [2 4 3 5] means ICA model 1 has 4 sources with
%             2,4,3,5 components respectively
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
% tol         Convergence tolerance     (Default = 0.0001)   
% max_steps   max no_steps to cease     (Default = 500)
% isonoise    1: learn isotropic noise  (Default)
%             0: learn non-isotropic diagonal noise
% ARD:        1: Automatic Relevance Determination is ON   (Default)
%             0: Automatic Relevance Determination is OFF.
% POS_SRC     1: learn positive bases  (Default for 'e', 't')   
%             0: learn any bases       (Default for 'g')
% POS_MIX     1: learn positive mixing matrix (Default for 'e', 't') 
%             0: learn any mixing matrix      (Default for 'g')
% nets        If ICA models are known or initialised in advance, these can 
%             be passed over as a cell array
%
%
%
% -----------
% Output 
% -----------
%
% HMM is a data structure with the following fields:
%
% K            The number of HMM states
% obs          The cell array of ICA networks
% train        The training parameters
% priors       The priors parameters
% gammas       The posterior state probabilities
% P            The state transition matrix
% Pi           The intial state probabilities
% posts        The posterior parameters  
% state_path   The Viterbi path for the data
% energy       The negative free energy
% model           'vbHMM-ICA'
% source_type     Type of ICA source model
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
addpath ../utils ../vbMoICA ../vbICA2 ../MixModel1d ../Shared

global CHECK_PROGRESS;

%--------------------Error catching and default setting--------------------
if nargin<2
  error('Please pass all necessary parameters!');
end
  
if nargin<3
    obsmodel = 'g';
end
    
if nargin<4
  CHECK_PROGRESS = 0;
else
  CHECK_PROGRESS = SEE_PDFS;
end


if nargin<5
    tol = 1e-4;
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

if nargin<11
    isnets = 0;
else 
  isnets = 1;
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
  error('isonoise can only be ON (1) or OFF (0)!')
end

if sum(ARD == [0 1])~=1
  error('ARD can only be ON (1) or OFF (0)!')
end


%--------------------Error catching and default setting--------------------






%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if isnets
  hmm = hmminit(data,states,src_type,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets);
else
  hmm = hmminit(data,states,src_type,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX);
end
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::





%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
hmm = hmmtrain(data,hmm);
hmm.model = 'vbHMM-ICA'
hmm.source_type = source_type;
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::