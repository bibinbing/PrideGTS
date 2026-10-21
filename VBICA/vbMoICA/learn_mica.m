function [Fica,nets] = learn_mica(data,nets,eta,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX);
% [Fica,nets] = learn_mica(data,nets,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,eta);
% 
% Learning of the individual ICA networks in an ICA mixture model 
% using variational Bayes.
% 
% Called from 'vbmoica' and 'ICA_update' (vbHMM).
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        no_rows=no_sensors. no_columns=no_data points
% states      Cell array that contains number states per ICA source MoG 
%             eg STATES{1} = [2 4 3 5] means ICA model 1 has 4 sources with
%             2,4,3,5 components respectively
% eta         The component probabilities for each data vector
%
%
% Optional parameters: 
%
% tol         Convergence tolerance     (Default = 0.0001)   
% max_steps   max no_steps to cease     (Default = 500)
% isonoise    1: learn isotropic noise  (Default)
%             0: learn non-isotropic diagonal noise
% ARD:        1: Automatic Relevance Determination is ON   (Default)
%             0: Automatic Relevance Determination is OFF
% POS_SRC     1: learn positive bases  (Default for 'e', 't')   
%             0: learn any bases       (Default for 'g')
% POS_MIX     1: learn positive mixing matrix. (Default for 'e', 't') 
%             0: learn any mixing matrix.      (Default for 'g')
%
%
% -----------
% Output 
% -----------
%
% MIXNET is a data structure with the following fields:
%
% energy        The free energy for the model 
% nets          The cell array containing the component ICA networks 
%               See 'help vbica2' for data structure
% priors.rho_0  The mixing prior hyperparameter before learning
% posts.rho     The mixing prior hyperparameter after learning
% eta           The component probabilities for each data vector
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
global CHECK_PROGRESS;
%--------------------Error catching and default setting--------------------
if nargin<3
  error('Please pass all necessary parameters!')
end

if nargin<4
  tol = 1e-5;
end

if nargin<5
  max_steps = 500;
end

if nargin<6
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



if tol<eps
  error('Please choose positive TOL!')
end

if max_steps<1
  error('Please choose MAX_STEPS > 1!')
end

if mod(max_steps,1)>0
  error('Please choose integer MAX_STEPS!')
end
 
if sum(ARD == [0 1])~=1
  error('ARD can only be ON (1) or OFF (0)!')
end
%--------------------Error catching and default setting--------------------




comps = length(nets);

for i=1:comps
  fprintf('Network %i\n',i)
  
  if CHECK_PROGRESS
    figure(i);
  end
  
  nets{i} = learn_ica2(data,nets{i},tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX,eta(i,:));

  Fica(i) = nets{i}.energy(length(nets{i}.energy));

end


