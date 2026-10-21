function [x]=sampgauss(m,C,N)
%
%  x=SAMPGAUSS(m,C,N)
%
%  samples N-times from an multi-dimensional gaussian distribution 
%  with covariance matrix C and mean m. Dimensionality is implied
%  in the mean vector
%
%  e.g: C=[1 .7;0.7 1];
%       m=[0;0];
%       x=sampgauss(m,C,300);
%(see e.g. B.D. Ripley, Stochastic Simulation, Wiley, 1987, pp. 98--99)
%
m=m(:);

r=size(C,1);
if size(C,2)~= r
  error('Wrong specification calling normal')
end
% find cholesky decomposition of A
[L,p]=chol(C);
% generate r  independent N(0,1) random numbers
z=randn(r,N);
x=m(:,ones(N,1))+L*z;

