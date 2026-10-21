function mixnet = vbmoica(data,states,source_type,SEE_PDFs,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,nets)
% mixnet = vbmoica(data,states,source_type,SEE_PDFs,tol,max_steps,ARD,POS_SRC,POS_MIX,nets)
% 
% ICA mixture model trained by Variational Bayes inc. 
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
% MIXNET is a data structure with the following fields:
%
% nets          The cell array containing the component ICA networks 
%               See 'help vbica2' for data structure
% priors.rho_0  The mixing prior hyperparameter before learning
% posts.rho     The mixing prior hyperparameter after learning
% eta           The component probabilities for each data vector
% energy_path   The negative free energy per iteration
% energy        The final NFE for the model 
% model         'vbMoICA'
% source_type   Type of ICA source model
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)


global CHECK_PROGRESS;

%--------------------Error catching and default setting--------------------
if nargin<2
  error('Please pass all necessary parameters!')
end

if nargin<3
  source_type = 'g';
end

if length(source_type)>1
  src_type = source_type(2);
else
  src_type = source_type;
end

if nargin<4
  CHECK_PROGRESS = 0;
else
  CHECK_PROGRESS= SEE_PDFs;
end


if nargin<5
  tol = 1e-5;
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
  INITIALISE = 1;
else
  INITIALISE = 0;
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


comps = length(states);
[sensors,points] = size(data);
% HARD_LIMIT is strict partitioning
HARD_LIMIT = 0;



%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
fprintf('Initialising...\n')
F = 0;
rho_0 = 5*ones(comps,1);



%If ICA models are not known, then initialise these as well
if INITIALISE
    [eta,nets] = initialise_mica(data,states,source_type,isonoise);
    if HARD_LIMIT
	eta = (eta==repmat(max(eta),comps,1))+eps;
    end
else
    eta = update_eta(data,nets,rho_0);
end

%plot probabilities
plot(eta')
axis([0 points -0.2 1.2])
drawnow



%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                INITIALISATION
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::











%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
err = tol+1;
steps = 0;
while err>tol
  steps = steps + 1;
  
  
  
  
  %---------------------------------M_Step---------------------------------
  %Update prior mix coefficients
  rho = rho_0 + sum(eta,2);
  %Update component ICA models
  [Fica nets] = learn_mica(data,nets,eta,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX);
  %---------------------------------M_Step---------------------------------
  
  
  

  
  
  %---------------------------------Energy---------------------------------
  %Work out free energy
  ent_eta = -sum(sum(eta.*log(eta+eps)));
  kl_bit1 = sum(gammaln(rho) - gammaln(rho_0));
  kl_bit2 = gammaln(sum(rho))-gammaln(sum(rho_0));
  kl_bit = kl_bit1-kl_bit2;
  old_F = F(length(F));
  F(steps) = sum(Fica)+kl_bit+ent_eta;
  %---------------------------------Energy---------------------------------
  
  
  
  
  
  %---------------------------------E_Step---------------------------------
  fprintf('\nUpdating mixture coefficients...\n')
  oldeta = eta(:);
  eta = update_eta(data,nets,rho);
  %erreta = norm(oldeta-eta(:));
  err = abs((F(steps)-old_F)/F(steps));
  %---------------------------------E_Step---------------------------------
  
  
  
  
  fprintf('\n[%i/%i]..........F_tot=%1.4f\n',steps,max_steps,F(steps));

  %plot probabilities
  plot(eta')
  axis([0 points -0.2 1.2])
  drawnow
  
  
  
  if steps == max_steps
    fprintf('\n\nReached maximum of %i training steps.\n',max_steps);
    break;
  end
end
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::





if steps<max_steps
  pcent = tol*100;
  fprintf('\n\nMixings converged to within %f percent in %i steps.\n',pcent,steps);
end





mixnet.nets = nets;
mixnet.priors.rho_0 = rho_0;
mixnet.posts.rho = rho;
mixnet.kappa = rho./sum(rho);
mixnet.eta = eta;

mixnet.energy_path = F;
mixnet.energy = F(length(F));

mixnet.model = 'vbMoICA';
mixnet.source_type = source_type;
