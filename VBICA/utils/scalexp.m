function eta = scalexp(log_eta,ALL);
% eta = scalexp(log_eta,ALL)
%
% Exponentiates log_eta while avoiding machine precision.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% log_eta   M x T data matrix
%
%
% Optional parameters
%
% ALL       0: Scales wrt largest element in row   (Default)
%           1: Scales wrtlargest element in matrix 
%
% -----------
% Output
% -----------
% 
% eta       Exponentiated log_eta
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 


if nargin <2
  ALL=0;
end
  [comps points] = size(log_eta);

if ALL
  scale = max(max(log_eta));
  index = log_eta-repmat(scale,comps,comps);
  bit1 = sum(sum(exp(index)));
  z = log(bit1)+scale;
  z = repmat(z,comps,comps);
  eta = exp(log_eta-z);
else
  scale = max(log_eta);
  index = log_eta-repmat(scale,comps,1);
  bit1 = sum(exp(index),1);
  z = log(bit1+eps)+scale;
  z = repmat(z,comps,1);
  eta = exp(log_eta-z);
end
