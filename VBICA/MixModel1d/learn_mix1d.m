function src = learn_mix1d(src,x,x_sq,tol,max_steps)
% src = learn_mix1d(src,x,x_sq,tol,max_steps)
%
% Train a 1-dimensional mixture model using the 
% Variational Bayes framework.
%
% Called from 'mixmodel1d', 'learn_ica1' and 'learn_ica2'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% src       Current model
% x         mixmodel1d: The observation data (sent as mxN)
%           vbICA1: Expectation of source signals (mxN) 
%           vbICA2: Location parameters of source posterior
% x_sq      mixmodel1d: x.^2 (sent as mxN)
%           vbICA1: Expectation of (source signal.^2) (mxN)
%           vbICA2: Precisions of source posterior
%
%
% Optional parameters
%
% tol       Convergence tolerance           (Default = 1e-5)
% max_steps Max number of iteration steps   (Default = 500)
%
%
% Note:     Exponentials take longer to converge, particularly 
%           if source pdfs not well described by mixture of   
%           exponentials. Therefore, tol=1e-6 for MoE. Truncated 
%           Gaussians are more robust and have better convergence.
%
%
%
% -----------
% Output
% -----------
%
% SRC is a data structure with the following fields:
%
% type             'g', 'e', 't' - Mixture of Gaussians, 
%                  exponentials or truncated Gaussians
%                  'hg','he','ht' - HMM versions
% m                The number of comp
% f_hidd           The negative free energy of hidden vars.
% f_params         The negative free energy of the params.
% fm               The total negative free energy
%
%                  In the field priors:
% lambda_0         Dirichlet parameters for mixing coeffs
% b_0,c_0          Gamma parameters for precisions
% m_0,tau_0        Normal parameters for means
%
%                  In the field posts:
% lambda           Dirichlet parameters  for mixing coeffs
% b,c              Gamma parameters for precisions
% mm,tau           Normal parameters for means
%
%                  Expected posterior values:
% pi               Mixing coefficients
% centres          Means 
% precs            Precisions
% gammas           Component probabilities
%
%
% If HMM, also
%
% pi               Initial state probabilities.
% P                Calculated trans. probs. for HMM
% posts.eta        Dirichlet parameters for P
% posts.lambda     Dirichlet parameters for pi
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)



global CHECK_PROGRESS;


if nargin<4
  tol = 1e-5;
end

if nargin<5
  max_steps = 500;
end

src_type = src.type;
pdf_fact = 0.5;
HMM=0;
if length(src_type) == 2
  HMM=1;
  src_type = src_type(2);
end




  


%========================Extract appropriate variables=====================
[m N]=size(x);  %N is number of data points
if m==1
  ALGORITHM = 1;
  m=src.m;        %m is the number of components in mixture model
  x = repmat(x,m,1);
  x_sq = repmat(x_sq,m,1);
  %Dummy variables for gammas.m
  m_q = x;
  b_q = x_sq;
elseif m==src.m
  ALGORITHM = 2;
  %Dummy variables for gammas.m
  m_q = x;
  b_q = x_sq;
  %Calculate 'true' x,x_sq.
  switch src_type
   case 'g'
    [x,x_sq] = deal(m_q,m_q.^2+1./b_q);
   otherwise
    [x,x_sq] = rect_expect(m_q,b_q);
  end
else
  error('Source signal dimensionality and source model mismatch!');
end
if(src_type == 'e')
  x_sq = x;
  pdf_fact=1;
  tol = min(tol,1e-6);
end

% PRIORS
% mixture prior
lambda_0=src.priors.lambda_0;
% component precision priors
b_0=src.priors.b_0;
c_0=src.priors.c_0;
% component mean priors
m_0=src.priors.m_0;
tau_0=src.priors.tau_0;



% POSTERIORS
% component precision posts
b=src.posts.b;
c=src.posts.c;
% component mean posts
mm=src.posts.mm;
tau=src.posts.tau;

%========================Extract appropriate variables=====================




% Initialise
ftot=0;
lik=[];
littlebit=eps;
Fgauss = 0;

% ALGORITHM = 2 for vbICA2. E-step only once for good convergence.
outer_steps = max_steps;
inner_steps = 1;
if ALGORITHM == 2
  outer_steps = 1;
  inner_steps = max_steps;
end

% Run...
for steps1=1:outer_steps



  %==================================E_Step================================
  if m==1
    gamma = ones(1,N);
  else
    [gamma xi] = gammas(src,m_q,b_q,ALGORITHM);
  end
  %==================================E_Step================================
  
  
  
  
  
  for steps2=1:inner_steps

  %==================================M_Step================================
  
  
  
  %--------------------Update lambda etc.--------------------
  if HMM
    
    % M-step for HMM specific variables 
    [src,Fp,Fpi,ent_gam] = hmm_mstep(src,gamma,xi);
    gamma_sum = sum(gamma');
    lambda = gamma_sum;  % for plotting only
    
    % contribution to energy    
    Fdir = Fp+Fpi;
    
  else
    
    gamma_sum = sum(gamma');
    lambda = lambda_0+gamma_sum;
    src.pi = lambda./sum(lambda);
    
    % store for E-step
    src.posts.lambda=lambda;
    
    % contribution to energy    
    lambda_p=lambda_0*ones(1,m);
    dir1 = sum(gammaln(lambda+eps) - gammaln(lambda_p+eps));
    dir2 = gammaln(sum(lambda+eps)) - gammaln(sum(lambda_p+eps));
    Fdir = dir1-dir2;
    ent_gam=-sum(sum(gamma.*log(gamma+eps)));
  end
  %--------------------Update lambda etc.--------------------
  
  
  
  %--------------------Update precisions---------------------
  mean_xsq = sum(gamma.*x_sq,2);
  mean_x = sum(gamma.*x,2);
  mu_sq = (gamma_sum.*(mm.^2+1./tau))';
  
  data_bit = mean_xsq-2*mm'.*mean_x+mu_sq;
  b = 1./(1/b_0+pdf_fact*data_bit');
  c = c_0+pdf_fact*gamma_sum;
  mean_beta = b.*c;
  
  % store for E-step
  src.posts.b=b;
  src.posts.c=c;
  
  % contribution to energy
  beta1 = gammaln(c) - gammaln(c_0);
  beta2 = c.*log(b) - c_0.*log(b_0);
  Fbeta = sum(beta1 + beta2);
  %--------------------Update precisions---------------------
  
  
  
  
  %-----------------------Update means-----------------------
  if src_type == 'g'
    tau = tau_0+mean_beta.*gamma_sum;
    mm = 1./tau.*(m_0+mean_beta.*mean_x');
    
    % store for E-step
    src.posts.mm=mm;
    src.posts.tau=tau;
    
    % contribution to energy
    b_ratio = tau_0./tau;
    gauss1 = -log(b_ratio);
    gauss2 = b_ratio-1;
    gauss3 = tau_0.*(mm-m_0).^2;
    Fgauss = -0.5*sum(gauss1 + gauss2 + gauss3);
  end
  %-----------------------Update means-----------------------
  

  
  %==================================M_Step================================


  

 
  
  
  
  
  %==================================Energy================================
  oldfm = ftot;
  f_hidd = ent_gam-N/2*log(2*pi);
  f_params =  Fgauss+Fbeta+Fdir;
  ftot = f_hidd+f_params;
  %==================================Energy================================

 

  % Plot if required
  if CHECK_PROGRESS
    switch src_type
     case 'g'
      plotMoG(mm,1./(b.*c),lambda./sum(lambda))
      drawnow
     case 'e'
      plotMoE(1./(b.*c),lambda./sum(lambda),max(max(x'))+1/max(sqrt(b.*c)))
      drawnow
     case 't'
      plotRMoG(1./(b.*c),lambda./sum(lambda),max(max(x'))+1/max(sqrt(b.*c)))
      drawnow
    end
  end
  
  
  % Convergence criterion for inner_steps
  err = abs((oldfm-ftot)/ftot);
  if err<tol
    break
  end
  
  end % End of inner_steps
  
  % Convergence criterion for outer_steps
  if err<tol
    break
  end

end  % End of outer_steps
  
% Put variables into data structure

src.centres=mm;
src.precs=mean_beta;
src.gammas=gamma;
src.f_hidd = f_hidd;
src.f_params = f_params;
src.ftot=ftot;
  




