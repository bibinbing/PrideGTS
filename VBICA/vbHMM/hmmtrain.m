function hmm=hmmtrain(x,hmm)
% hmm=hmmtrain(x,hmm)
%
% Train vbHMM-ICA model using variational Baum Welch algorithm.
%
% Called from 'vbhmm'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% x         No. rows = No. sensors, No. columns = No. data points
% hmm       Initialised vbHMM-ICA model
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
% 
%
% --------------------------------------------------------------
%
% Original code by Iead Rezek
% Paper: Learning Ensemble Hidden Markov Models for Biosignal 
%        Analysis (www.robots.ox.ac.uk/~parg)
%
% Modified by Rizwan Choudrey for use in vbICA model
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)



x = x';


p=length(x(1,:));
T=length(x(:,1));


if isfield(hmm,'K')
  K=hmm.K;
else
  disp('Error in hmmtrain: K not specified');
  return
end

if ~isfield(hmm,'train')
  disp('Error in hmmtrain: hmm.train not specified');
  return
end

if isfield(hmm.train,'cyc')
  cyc=hmm.train.cyc;
else
  cyc=100; 
end

if isfield(hmm.train,'tol')
  tol=hmm.train.tol;
else
  tol=0.0001; 
end

if ~isfield(hmm.train,'obsupdate')
  updateobs=ones(1,hmm.K);  % update observation models for all states
else
  updateobs=hmm.train.obsupdate;
end

if ~isfield(hmm.train,'pupdate')
  updatep=1;
else
  updatep=hmm.train.pupdate;
end


if ~isfield(hmm,'pi')
  Pi=rand(K,1);
  Pi=Pi/sum(Pi);
else
  lambda_hat = hmm.posts.lambda_hat;
  sumlambda = sum(lambda_hat);
  for i=1:K
      index = digamma(lambd_hat(i))-digamma(sumlambda);
      Pi(i,1) = exp(index);
  end
end


if ~isfield(hmm,'P')
  P=rand(K);
  P=rdiv(P,rsum(P));
else
  eta_hat = hmm.posts.eta_hat;
  etasum = sum(eta_hat,2);
  for i=1:K
      for j=1:K
	  index = digamma(eta_hat(i,j))-digamma(etasum(i));
	  P(i,j) = exp(index);
      end
  end
end

LL=[];
lik=0;

alpha=zeros(T,K);
beta=zeros(T,K);
if ~isfield(hmm,'gammas')
  Gamma = rand(T,K);
  Gamma = Gamma./repmat(sum(Gamma,2),1,K);
  Gammasum = sum(Gamma);
else
  Gamma = hmm.gammas;
  Gammasum = sum(Gamma);
end

  gamma=zeros(T,K);
F=0;

for cycle=1:cyc
  
  %%%% FORWARD-BACKWARD 
   % Observation model

   if sum(updateobs) > 0
     hmm = ICA_update(x',Gamma',hmm);
   end
  
  
  
  
  Gamma=[];
  Gammasum=zeros(1,K);
  Scale=zeros(T,1);
  Xi=zeros(T-1,K*K);
  
  for n=1:1
    
    B = ICA_like(x',hmm);
    scale=zeros(T,1);
    alpha(1,:)=Pi'.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for i=2:T
      alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
      scale(i)=sum(alpha(i,:));
      alpha(i,:)=alpha(i,:)/scale(i);
    end;
    
    beta(T,:)=ones(1,K)/scale(T);
    for i=T-1:-1:1
      beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i); 
    end;
    
    oldgamma = gamma;
    
    gamma=(alpha.*beta);
    gamma=rdiv(gamma,rsum(gamma));
    gammasum=sum(gamma);
    err = norm(gamma-oldgamma);
    
    xi=zeros(T-1,K*K);
    for i=1:T-1
      t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
      xi(i,:)=t(:)'/sum(t(:));
    end;
    
    Scale=Scale+log(scale);
    Gamma=[Gamma; gamma];
    Gammasum=Gammasum+gammasum;
    Xi=Xi+xi;
  end;
  
  

  
  %%%% M STEP 

  % transition matrix 
  sxi=sum(Xi);
  sxi=reshape(sxi,K,K);
  eta0 = hmm.priors.eta_0;
  eta_hat = eta0+sxi;
  seta = repmat(sum(eta_hat,2),1,K);
  if updatep
    trueP=eta_hat./seta;
  end
  for i=1:K
      for j=1:K
	  index = digamma(eta_hat(i,j))-digamma(seta(i));
	  P(i,j) = exp(index);
      end
  end
  
  %contribution to free energy
  if length(eta0(1,:))~=K
      eta0 = repmat(eta0,1,K);
  end
  bit1 = sum(gammaln(eta_hat) - gammaln(eta0),2);
  bit2 = gammaln(sum(eta_hat,2))-gammaln(sum(eta0,2));
  Fp = sum(bit1-bit2);
  
  
  
  
  % initial state probabilities
  lambda0 = hmm.priors.lambda_0;
  lambda_hat = lambda0+Gamma(1,:)';
  sumlambda = sum(lambda_hat);
  truePi = lambda_hat./sum(lambda_hat);
  for i=1:K
      index = digamma(lambda_hat(i))-digamma(sumlambda);
      Pi(i,1) = exp(index);
  end
   
  %contribution to free energy
  bit1 = sum(gammaln(lambda_hat) - gammaln(lambda0));
  bit2 = gammaln(sum(lambda_hat))-gammaln(sum(lambda0));
  Fpi = sum(bit1-bit2);

  Fstate = Fp+Fpi;

 

 
  
  entStart = -sum(Gamma(1,:).*log(Gamma(1,:)+eps));
  for i=1:K
    Xi3d = reshape(Xi,[T-1,K,K]);
    sxi = sum(Xi3d(:,:,i),2);
    Psi(:,:,i) = Xi3d(:,:,i)./(repmat(sxi,1,K)+eps);
  end
  entXi = -sum(sum(sum(Xi3d.*log(Psi+eps))));
  entHidden = entStart+entXi;
  

  oldF=F;
  F = Fstate+hmm.obs.energy+entHidden;
  err = abs((F-oldF)/F);

  fprintf('cycle %i Ftot = %f ',cycle,F);  
  
  if err<tol
  fprintf('\n');
    break;
  end;
  fprintf('\n');

end



hmm.P=trueP;
hmm.Pi=truePi;

hmm.posts.lambda_hat = lambda_hat;
hmm.posts.eta_hat = eta_hat;

hmm.gammas = gamma';

[block,ll,llb] = hmmdecode(x',hmm);
hmm.state_path = block.q_star;
hmm.energy = F;


