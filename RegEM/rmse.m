% RMSE - Evaluates root mean squared error.
%
% E = RMSE(X)
% 
% E = RMSE(X1,X2)

function e = rmse(X, Xh)

if nargin < 2
  e = sqrt(X(:)'*X(:)/numel(X));
else
  z = X(:)-Xh(:);
  e = sqrt(z'*z/numel(z));
end
