function log_like = log_ptilde1(y,y_sq,src);
% log_like = log_ptilde1(y,y_sq,src)
%
% Calculates the pseudo log-likelihood of y under src. Used
% in 'mixmodel1d' and 'vbica1'.
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
% y         mixmodel1d: The observation data (sent as mxN)
%           vbICA: Expectation of source signals (mxN) 
% y_sq      mixmodel1d: y.^2 (sent as mxN)
%           vbICA: Expectation of (source signal.^2) (mxN)
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

b=src.posts.b;
c=src.posts.c;
mean_beta=(b.*c)';
if m~=length(mean_beta)
  m = length(mean_beta);
  y = repmat(y,m,1);
  y_sq = repmat(y_sq,m,1);
end

src_type = src.type;
if length(src_type) == 2
  src_type = src_type(2);
end

pdf_fact = 0.5;
mm_sq = zeros(m,1);
tgauss_fact = zeros(m,N);

%Adjust for source type
switch src_type
 case 'g'
  data_term = y_sq;
  mm = src.centres';
  tau = src.posts.tau';
  mm_sq = mm.^2+1./tau;
 case 't'
  data_term = y_sq;
  mm = mm_sq;
 case 'e'
  data_term = y;
  pdf_fact = 1;
  mm=mm_sq;
end



for s=1:m,
   log_tilde_beta(s)=digamma(c(s))+log(b(s)+eps);
end

data_bit = -pdf_fact*repmat(mean_beta,1,N).*(data_term+repmat(mm_sq, ...
						  1,N)-2*repmat(mm,1,N).*y);
log_like = pdf_fact*repmat(log_tilde_beta',1,N) + data_bit;

