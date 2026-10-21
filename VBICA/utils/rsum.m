function Z=rsum(X)
% Z=rsum(X)
%
% row sum
%
%
% ------------------------------------------------------
%
% Original code by Iead Rezek

  Z=zeros(size(X(:,1)));

for i=1:length(X(1,:))
  Z=Z+X(:,i);
end
