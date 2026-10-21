function VBRPCAResult = RunVBRPCA(X,NumComponents,NumReconComp,Maxiters,verbose,initial)

VBRPCAResult = struct;
options = struct;
options.MAXITER = Maxiters;
options.init = initial;
if isempty(NumComponents)
    options.initial_rank = 'auto';
else
    options.initial_rank = NumComponents;
end
options.verbose = verbose;

[X_hat, A_hat, B_hat, E_hat, gamma,alpha,beta,it,t,Sigma_A,Sigma_B] = VBRPCA(X,options);

sigma2 = 1 ./ gamma; % gamma is the prior variance of principle components.
total_sigma2 = sum(sigma2);
Contribution = (sigma2 / total_sigma2).*100;
[SortContri, Idx] = sort(Contribution, 'descend');
A_hat = A_hat(:,Idx);
B_hat = B_hat(:,Idx);

Components = A_hat;
SpatialResponses = B_hat;
cme = A_hat(:,1:NumReconComp) * B_hat(:,1:NumReconComp)';
FilteredData = X - cme;

VBRPCAResult.Contribution = SortContri;
VBRPCAResult.FilteredData = FilteredData;
VBRPCAResult.OriginTimeSeries = X;
VBRPCAResult.SpatialResponsesMatrix = SpatialResponses; % ATTENTION!!!
VBRPCAResult.Components = Components; % Store sorted components

end

