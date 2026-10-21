function log_like = log_ptilde2(y,y_prec,src);
% log_like = log_ptilde2(y,y_prec,src)
%
% Calculates the pseudo log-likelihood of y under src. Used
% in 'vbica2'.
%
% Called from 'gammas.m', 'hmmdecode_ica', 'FB'. 
%
%
% -----------
% Input
% -----------
%
% Necessary parameters
%
% y         Location parameters of source signal
% y_prec    Precisions of source signal
% src       Current model
%
%
%
% -----------
% Output
% -----------
%
% log_like  Calculated pseudo log-likelihood
%
%
% --------------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
[m,N] = size(y);
littlebit = eps;

b = src.posts.b;
c = src.posts.c;

src_type = src.type;
if length(src_type) == 2
  src_type = src_type(2);
end

pdf_fact = 0.5;
mm_sq = zeros(m,1);
tgauss_fact = zeros(m,N);

switch src_type
 case 'g'
  mm = src.centres;
  tau = src.posts.tau;
  mm_sq = (mm.^2+1./tau)';
  yfact = y;
 case 't'
  yfact = y;
  tgauss_fact = log(erfc(-0.5*y.*y_prec)+eps);
 case 'e'
  yfact = 0.5*rect_expect(y,y_prec);
  pdf_fact=1;
end




for s=1:m
   log_tilde_beta(s)=psi(c(s))+log(b(s)+eps);
end

log_beta_tilde = repmat(log_tilde_beta',1,N);
mean_beta=c.*b;
databit = pdf_fact*(y.*yfact.*y_prec - repmat((mm_sq.*mean_beta'),1,N) - ...
		    log(y_prec));
log_like = pdf_fact*log_beta_tilde + databit+tgauss_fact;
