function log_like = avge_like(data,net,x,x_prec);
% log_like = avge_like(data,net,x,x_prec);
%
% Finds average log likelihood of data under given vbICA model.
% If source posterior parameters not given explicitely, then uses 
% the ones stored in net.
%
% Called from 'update_eta' (vbmoica) and 'ICA_like' (vbhmm).
%
%
% -----------
% Input
% -----------
%
% Necessary parameters:
%
% data        No. rows = No. sensors, No. columns = No. data points
% network     ICA network
%
%
% Optional parameters
%
% x           Expectation of source signals
% x_prec      Precision of expectations
%
%
%
% -----------
% Output
% -----------
%
% log_like    Average log likelihood of data under net
%
%
% ------------------------------------------------------
%
% Original code by Rizwan Choudrey 
% Thesis: Variational Methods for Bayesian Independent
%         Component Analysis (www.robots.ox.ac.uk/~parg)
  
if nargin<4
    x = net.recon;
    x_prec = net.posts.x_prec;
end

    
H = net.mixmatrix;
H_prec = net.posts.mix.H_prec;


noise = net.noise;
mn = net.posts.noise.mn;
mnb = net.posts.noise.mb;
c_Lam_hat = mean(net.posts.noise.c);
b_Lam_hat = mean(net.posts.noise.b);

    
[sensors points] = size(data);
[sensors sources] = size(H);
dc = repmat(mn,1,points); 

switch net.source_type
 case 'g'
  H_sq = H.^2+1./H_prec;
 otherwise
  H_sq = net.posts.mix.H_mean.*H+1./H_prec;
end

x_sq = x.^2+1./x_prec;
  
data_hat = H*x;
dd=(data-dc).^2;
ddh=(data-dc).*data_hat;

data_hat_var = H_sq*x_sq-(H.^2)*(x.^2);
dhdh = data_hat.^2+data_hat_var;

mean_bit = repmat(1./mnb,1,points);
Ed = dd-2*ddh+dhdh+mean_bit;

prelike = repmat(-0.5*noise,1,points).*Ed;
for i=1:length(b_Lam_hat)
  log_lam_tild(i) = psi(c_Lam_hat(i))+log(b_Lam_hat(i));
end
if length(log_lam_tild)<sensors
  log_lam_tild = repmat(log_lam_tild,1,sensors);
end

log_like = sum(repmat(0.5*log_lam_tild',1,points) + prelike,1);

    
    
    
    
    
    
    
