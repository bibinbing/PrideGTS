function hmm = hmminit (x,states,source_type,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets)
% hmm = hmminit (x,states,source_type,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets)
%
% Initialise  vbHMM-ICA model. 
%
% Called from 'vbhmm'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% x           No. rows = No. sensors, No. columns = No. data points
% states      Cell array that contains number states per ICA source MoG 
%             eg STATES{1} = [2 4 3 5] means ICA model 1 has 4 sources with
%             2,4,3,5 components respectively.
%
%
% Optional parameters: 
%
% source_type Mixture of 'g'aussians     (Default)
%             Mixture of 'e'xponentials
%             Mixture of 't'runcated gaussians
%             'hg','he','ht' - HMM versions
% tol         Convergence tolerance     (Default = 0.0001)   
% max_steps   Max no_steps to cease     (Default = 500)
% isonoise    1: learn isotropic noise  (Default)
%             0: learn non-isotropic diagonal noise
% ARD:        1: Automatic Relevance Determination is ON   (Default)
%             0: Automatic Relevance Determination is OFF
% POS_SRC     1: learn positive bases  (Default for 'e', 't')   
%             0: learn any bases       (Default for 'g')
% POS_MIX     1: learn positive mixing matrix (Default for 'e', 't') 
%             0: learn any mixing matrix.     (Default for 'g')
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
% K      The number of HMM states
% obs    The cell array of intialised ICA networks
% train  The training parameters
% priors The prior parameters
% gammas The initialised posterior state probabilities
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
if nargin< 3
  source_type = 'g';
end

if length(source_type)>1
  src_type = source_type(2);
else
  src_type = source_type;
end

if nargin < 4
  tol = 1e-4;
end

if nargin < 5
  max_steps = 500;
end

if nargin < 6
  isonoise = 1;
end

if nargin<7
  ARD = 1;
end


if nargin<8
  switch src_type
   case 'g'
    POS_SRC = 0;
   otherwise
    POS_SRC = 1;
  end
end

if nargin<9
  switch src_type
   case 'g'
    POS_MIX = 0;
   otherwise
    POS_MIX = 1;
  end
end

if nargin<10
  INITIALISE = 1;
else
  INITIALISE = 0;
end

T = length(x);

if iscell(states)
    hmm.K = length(states);
else
    hmm.K = states;
    states{1} = states;
end
hmm.obs.states = states;
hmm.obs.isonoise = isonoise;  

hmm.obs.ARD = ARD;
hmm.obs.POS_SRC = POS_SRC;
hmm.obs.POS_MIX = POS_MIX;


hmm.train = 1;
hmm.train.tol = tol;
hmm.train.cyc = max_steps;


pseudo_count_lambda = 1;
hmm.priors.lambda_0 = repmat(pseudo_count_lambda,hmm.K,1);
base_count = 5;
same_comp_count = 5;
hmm.priors.eta_0 = base_count*ones(hmm.K);
hmm.priors.eta_0 = hmm.priors.eta_0+same_comp_count*eye(hmm.K);


if INTIALISE
  [eta nets] = initialise_mica(x,states,source_type,isonoise);
else
  eta = update_eta(x,nets,hmm.priors.lambda_0);
end

hmm.gammas=eta';
hmm.obs.ica = nets;


