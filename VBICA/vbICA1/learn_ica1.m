function network = learn_ica1(data,network,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX)
% network = learn_ica1(data,network,tol,max_steps,isonoise,ARD,POS_SRC,POS_MIX)
%
% Learn ICA model with Automatic Relevance Determination (ARD) usng 
% Variational Bayes with *Gaussian* source posterior.
% 
% y = Hx + e
%
% where each x has a MoG prior density and e is Gaussian noise.
% ARD scores each source on relevance and suppresses unnecessary sources.
%
% Called from 'vbica1'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        No. rows = No. sensors, No. columns = No. data points
% network     Initialised ICA network ready for training
%
%
% Optional parameters
%
% tol         Convergence tolerance    (Default = 1e-4)
% max_steps   Max no_steps to cease    (Default = 500)
% isonoise    1: learn isotropic noise (Default)
%             0: learn non-isotropic diagonal noise
% ARD         1: Automatic Relevance Determination is ON  (Default)
%             0: Automatic Relevance Determination is OFF
% POS_SRC     1: learn positive bases  (Default for 'e', 't')   
%             0: learn any bases.      (Default for 'g')
% POS_MIX     1: learn positive mixing matrix. (Default for 'e', 't') 
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
  
  
  
%-------------------Error catching and default setting---------------------
if nargin<2
  error('Please pass all necessary parameters!')
end

if nargin<3
  tol = 0.0001;
end

if nargin<4
  max_steps = 500;
end

if nargin<5
  isonoise = 1;
end

if nargin<6
  ARD = 1;
end

source_type = network.source_type;   
if length(source_type)==2
  src_type = source_type(2);
else
  src_type = source_type;
end

if nargin<7
  switch src_type
   case 'g'
    POS_SRC = 0;
   otherwise
    POS_SRC = 1;
  end
end

if nargin<8
  switch src_type
   case 'g'
    POS_MIX = 0;
   otherwise
    POS_MIX = 1;
  end
end
%-------------------Error catching and default setting---------------------






%----------------------Extract appropriate variables-----------------------
global CHECK_PROGRESS;

source = network.sources;             % Cell of source models
states = network.source_states;       % Vector of states per source
priors = network.priors;              % Structure holding priors
alpha = network.alphas;               % ARD coefficients
mu_x_hat = network.posts.x_mean;      % Initialised source posterior means
beta_x_hat = network.posts.x_prec;    % Initialised source posterior precs
m_H_hat = network.posts.mix.H_mean;   % Initialised mixmatrix post. means
alpha_hat = network.posts.mix.H_prec; % Initialised mixmatrix post. precs
lambda_hat = network.noise;           % Initialised noise posterior precs
mn = network.posts.noise.mn;
mb = network.posts.noise.mb;

[sensors,points] = size(data);
sources = length(states);

MEAN_UPDATE=1;
f_dc=0;
dc = repmat(mn,1,points);


% Rectify accordingly if positivity priors
if POS_SRC
  MEAN_UPDATE=0;
  [x,x_sq] = rect_expect(mu_x_hat,beta_x_hat);
else
  [x,x_sq] = deal(mu_x_hat,mu_x_hat.^2+1./beta_x_hat);
end

if POS_MIX
  [H,H_sq] = rect_expect(m_H_hat,alpha_hat);
else
  [H,H_sq] = deal(m_H_hat,m_H_hat.^2+1./alpha_hat);
end
%----------------------Extract appropriate variables-----------------------





%flow control parameters
%(part of bug testing, but left in to allow modification if so wished)
factSource = 1;
fact2 = 1;
factx = 1;
facth = 1;
src_steps = max_steps;
if length(source_type)==2       % If HMM sources, otherwise REALLY slow!
  src_steps = min(5,max_steps); % However, NFE may not be monotonic 
end 
x_steps = max_steps;
h_steps = max_steps;
max2_steps = max_steps;
anneal=1;





%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


fprintf('Learning...\n\n')
steps=0;
err = tol+1;
ftot=0;



while err>tol
  steps=steps+1;

  fprintf('Updating sources...')
  %----------------------------source models-------------------------------
  src_tol = factSource*tol;
  for i=1:sources
    if CHECK_PROGRESS
      subplot(1,sources,i)
    end
    source{i} = learn_mix1d(source{i},x(i,:),x_sq(i,:),src_tol,src_steps);
    fsp(i) = source{i}.f_params;
    fsh(i) = source{i}.f_hidd+points/2*log(2*pi);
  end
  %----------------------------source models-------------------------------

  

 

    fprintf('Updating observation network...')
  
  %=====Loop estimation of non-MoG parameters to allow for convergence=====
  
  iter2 = 0;
  err2 = fact2*tol+1;
  
  while err2>fact2*tol
    iter2 = iter2+1;
    
    clear mu_x_hat beta_x_hat m_H_hat alpha_hat


    %--------------------------source signals------------------------------
    xtol = factx*tol;
    oldoldx = x(:);
    [x,x_sq,mu_x_hat,beta_x_hat,ent_x] = recon_source1(data-dc,x,H,H_sq,lambda_hat,source,xtol,x_steps,POS_SRC);
    errx2 = x(:)-oldoldx;
    %--------------------------source signals------------------------------

    
    switch src_type
     case 'e'
      alpha_prec = zeros(sources,1);
      alpha_mean = alpha;
     otherwise
      alpha_prec = alpha;
      alpha_mean = zeros(sources,1);
    end

    
    %----------------------------mixmatrix---------------------------------
    htol = facth*tol;
    oldoldh = H(:);
    [H,H_sq,m_H_hat,alpha_hat,f_h] = learn_matrix(data-dc,H,x,x_sq, ...
						  lambda_hat,alpha_prec,alpha_mean,htol,h_steps,POS_MIX);
    errh2 = H(:)-oldoldh;
    %----------------------------mixmatrix---------------------------------


    
    %------------------------ARD coefficients alpha------------------------
    % ARD doesn't work for 'e' or 'he', but left in for experimentation
    switch src_type
     case 'e'
      H_type = H;
     otherwise
      H_type = H_sq;
    end
    [b_alpha_hat,c_alpha_hat,f_alph] = update_alpha(priors.mix,H_type,src_type,ARD);
    alpha = b_alpha_hat.*c_alpha_hat;
    %------------------------ARD coefficients alpha------------------------

    
    if (iter2==max2_steps)
      break;
    end
    
    err2 = sqrt(sum(errx2.^2)+sum(errh2.^2));

  end  % End of network convergence loop
  
  %=====Loop estimation of non-MoG parameters to allow for convergence===== 


  
  

  
  if MEAN_UPDATE
  %--------------------------Sensor noise offset---------------------------
  [mn,mb,f_dc] = update_mean(priors.noise,data,H,x,lambda_hat);
  dc = repmat(mn,1,points);
  %--------------------------Sensor noise offset---------------------------
  end
  

  %-----------------------------Sensor noise-------------------------------
  [b_Lam_hat,c_Lam_hat,f_lam] = learn_noise(priors.noise,data-dc,mb,H,x,H_sq,x_sq,isonoise);
  lambda_hat = b_Lam_hat.*c_Lam_hat;
  %-----------------------------Sensor noise-------------------------------



  %-----------------------------Free Energy--------------------------------
  % Calculate free energy
  % Note: these formulae assume that the posteriors
  % have just been optimised simultaneously, therefore
  % ftot may not be monotonic when parameters are changing rapidly.
  
  f_mix = f_alph+f_h;
  f_noise = f_dc+f_lam;

  f_params = sum(fsp)+f_mix+f_noise;
  f_hidd = sum(fsh)+ent_x-0.5*sensors*points*log(2*pi);
  
  % Store old values
  oldf=ftot;
  
  % Put it together
  ftot = f_params+f_hidd;

  F(steps) = ftot;
  fprintf('\n[%i/%i]..........F_tot=%1.4f\n',steps,max_steps,ftot);
  err = abs((ftot-oldf)/ftot);
  %-----------------------------Free Energy--------------------------------

  tol = tol*anneal;
  
  if steps==max_steps
      fprintf('\n\nReached maximum of %i training steps.\n',max_steps);
    break;
  end
  
end  % End of overall loop


if steps<max_steps
  pcent = tol*100;
  fprintf('\n\nEnergy converged to within %f percent in %i steps.\n',pcent,steps);
end

%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%                                  LEARNING
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




fprintf('State paths...')
for i=1:sources
  if length(network.source_type)==2
    [block,ll,llb] = hmmdecode1(x(i,:),x_sq(i,:),source{i});
    state_path(i,:) = block.q_star; %HMM state paths
  else
    [p,state_path(i,:)] = max(source{i}.gammas);
  end
end
network.state_paths = state_path;
fprintf('\n')


network.recon = x;
network.sources = source;
network.mixmatrix = H;
network.alphas = alpha;
network.noise = lambda_hat;
network.mean = mn;

network.priors = priors;

network.posts.x_mean=mu_x_hat;
network.posts.x_prec=beta_x_hat;
network.posts.mix.H_mean=m_H_hat;
network.posts.mix.H_prec=alpha_hat;
network.posts.alpha.b=b_alpha_hat;
network.posts.alpha.c=c_alpha_hat;
network.posts.noise.b=b_Lam_hat;
network.posts.noise.c=c_Lam_hat;
network.posts.noise.mn = mn;
network.posts.noise.mb = mb;



network.energy_path = F;
network.f_params = f_params;
network.f_hidd = f_hidd;
network.energy = ftot;
network.model = 'vbICA1';

