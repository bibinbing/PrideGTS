function mix1 = kmeans1(k,y)
% mix1 = kmeans1(k,y)
%
% Train a 1-dimensional k-means model.
%
% Called from 'init_MoG'.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% k          Number of components
% y          Data 
%
%
%
% -----------
% Output
% -----------
%
% The fields in MIX1 are:
%
% k                The number of components
% m                Vector of means
% v                Vector of variances
% pi               Vector of mixing proportions
% nloops           Number of iterations used
%
%
%
% -------------------------------------------------------
%
% Original code by Will Penny
%
% Modified by Rizwan Choudrey for use in vbICA model
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

if nargin<3
  eta = 1;
end


y=y(:)';
N=length(y);
MIXTURES = 1;


if length(eta)==1
  MIXTURES = 0;
  eta = ones(1,N);
end


% Spread seeds evenly according to CDF
[x,i]=sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i=ones(1,N);
m=x(seeds);

for loops=1:100,  
  for j=1:k,
   d(j,:)=(y-m(j)).^2;
  end 

 [tmp,i]=min(d);
 if sum(i-last_i)==0
   % If assignment is unchanged
   break;
 else
   % Recompute centres
   for j=1:k,
      if sum(eta(i==j))==0,
	m(j)=sum(eta(i==j).*y(i==j));
      else
	m(j)=sum(eta(i==j).*y(i==j))/sum(eta(i==j));
      end;
   end
   last_i=i;
 end
end  

% Compute variances and mixing proportions
for j=1:k,
  v(j)=(sum(eta(i==j).*(y(i==j)-m(j)).^2))/(sum(eta(i==j))+eps);
  if v(j)==0
    v(j) = 1000;
  end
if MIXTURES
  mix_prob(j)=sum(eta(i==j));
else
  mix_prob(j)=length(y(i==j))/N;
end
  gammas(j,:) = 1/(2*pi*v(j)).*exp(-((y-m(j)).^2)./(2*v(j)))*mix_prob(j);
end
sumg = repmat(sum(gammas),k,1);
gammas = gammas./sumg;

mix1.gammas = gammas;
mix1.v=v;
mix1.m=m;
mix1.pi=mix_prob./sum(mix_prob);
mix1.k=k;
mix1.nloops=loops;
mix1.last_i=last_i;