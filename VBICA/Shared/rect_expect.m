function [x,x_sq] = rect_expect(m,b);
% [x,x_sq] = rect_expect(m,b);
%
% Works out expectations under truncated Gaussian density.
%
% Called from various functions.
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% m        Location parameter of density
% b        Precision of density
%
%
%
% -----------
% Output
% -----------
%
% x,x_sq  Expectations
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
k1 = sqrt(2./(pi*b));
k2 = 1./erfcx(-m.*sqrt(b/2));
pos_fact = k1.*k2;

x = m + pos_fact; 
if nargout > 1
  x_sq = m.*x + 1./b;
end
