% RMSE - Evaluates root mean squared error.
%
% E = RMSE(X)
% 
% E = RMSE(X1,X2)

function e = rmse(X, Xh)

if nargin < 2
  % RMS
   e = sqrt(sum(X(:).^2) / numel(X));
else
  % RMSE
  z = X(:)-Xh(:);
  e = sqrt(z'*z/numel(z));
end
