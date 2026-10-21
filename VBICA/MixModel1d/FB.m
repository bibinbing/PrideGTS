function [Gamma,Xi] = FB(Pi,B,P)
% [Gamma,Xi] = FB(Pi,B,P)
%
% This is the Forward-Backward step for an HMM. 
%
% Called from 'gammas' (MixModel1d). 
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% Pi        Initial state probabilities
% B         Observation likelihood
% P         Transition matrix
%
%
%
% -----------
% Output
% -----------
%
% Gamma     Calculated state probabilities
% Xi        Calculated transition probabilities
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


  % Initialise
  [T,K] = size(B);
  Gamma=[];
  Gammasum=zeros(1,K);
  Scale=zeros(T,1);
  Xi=zeros(T-1,K*K);

  % Forward step
  scale=zeros(T,1);
  alpha(1,:)=Pi.*B(1,:);
  scale(1)=sum(alpha(1,:));
  alpha(1,:)=alpha(1,:)/scale(1);
  for i=2:T
    alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
    scale(i)=sum(alpha(i,:));
    alpha(i,:)=alpha(i,:)/scale(i);
  end;
  
  % Backward step
  beta(T,:)=ones(1,K)/scale(T);
  for i=T-1:-1:1
    beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i); 
  end;
  
  % Calculate state probs.
  gamma=(alpha.*beta);
  gamma=rdiv(gamma,rsum(gamma));
  gammasum=sum(gamma);
  
  % Calculate trans probs.
  xi=zeros(T-1,K*K);
  for i=1:T-1
    t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
    xi(i,:)=t(:)'/sum(t(:));
  end;
  
  % Scale to avoid machine precision
  Scale=Scale+log(scale);
  Gamma=[Gamma; gamma];
  Gamma = Gamma';
  Xi=Xi+xi;

  