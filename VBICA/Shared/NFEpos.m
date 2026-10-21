function F = NFEpos(m,b);
% [Gamma,Xi] = FB(Pi,B,P)
%
% Calculates the negative free energy of a truncated
% Gaussian.
%
% Called from 'recon_source'. 
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
% F        Negative free energy
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
  
logerf = log(erfc(-m.*b./2)+eps);
erfbit = 1./erfcx( -m.*sqrt(b./2));
erfcoeff = m.*sqrt(b./(2*pi));
ent = -0.5*(log(b)+1);
F = ent-erfcoeff.*erfbit+logerf;