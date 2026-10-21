function y = scale01(x);
%   y = scale01(x);
%
% Scales data to lie within [0,1]
  
  y = x - min(min(x)) + eps;
  y = y/max(max(y));
  
return;
  
