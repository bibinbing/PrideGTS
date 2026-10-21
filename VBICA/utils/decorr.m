function [newdata,transform] = decorr(data,fact);
% [newdata,transform] = decorr(data,fact);
%
% Decorrelates data by varying amounts depending on fact.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% data      M x T data matrix
%
%
% Optional parameters
%
% fact      Power of decorrelation (Default = 1)
%
%
%
% -----------
% Output
% -----------
% 
% newdata   Decorrelated data
% transform Decorrelating matrix
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)

  
if nargin<2
  fact = 1;
end


[dim,points] = size(data);
ndata = normalis(data',data')';

for i=1:points
  cov3d(:,:,i) = ndata(:,i)*ndata(:,i)';
end
cov2d = sum(cov3d,3)./points;
[v,d] = eig(cov2d);
d = abs(d);
transform = v*inv(d).^(0.5*fact)*v';
newdata = transform*ndata;