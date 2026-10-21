function VBPCAResult = RunVBPCA(X, NumComponents, NumReconComp, HypUpData, MaxIters, Display)
%RUNVBPCA1 Perform Variational Bayesian PCA and process results.
%
% Inputs:
%   X - Input data matrix (n x p)
%   NumComponents - Number of principal components to extract
%   NumReconComp - Component index for reconstructing CME
%   HypUpData - Hyperparameter update data (e.g., number of iterations)
%   MaxIters - Maximum number of iterations for VBPCA
%   Display - Display intermediate results (true/false)
%
% Outputs:
%   VBPCAResult - Struct containing:
%       SpatialResponsesRatio - Normalized spatial responses (percent)
%       Contribution - Sorted contribution of each component
%       FilteredData - Data with CME removed
%       SpatialResponsesMatrix - Sorted components


% Initialize output struct
VBPCAResult = struct;

% Input validation
if isempty(X)
    error('Input data X cannot be empty.');
end
[n, p] = size(X);
if NumComponents > min(n, p)
    error('NumComponents cannot exceed the minimum dimension of X.');
end
if NumReconComp > NumComponents || NumReconComp < 1
    error('NumReconComp must be between 1 and NumComponents.');
end
if MaxIters < 1

    error('MaxIters must be a positive integer.');
end

% Run VBPCA
try
    Result = vbpca_full(X, NumComponents, 'display', Display, ...
        'niter_broadprior', HypUpData, 'maxiters', MaxIters);
catch err
    error('Error in vbpca_full: %s', err.message);
end

% Compute contribution and sort
Contribution = (Result.Va / sum(Result.Va)) * 100;
[SortContri, Idx] = sort(Contribution, 'descend');
Components = Result.A(:, Idx);
SpatialResponses = Result.S(Idx, :);

% Reconstruct CME and filter data
cme = Components(:, NumReconComp) * SpatialResponses(NumReconComp, :);
FilteredData = X - cme;

% Store results in struct
VBPCAResult.Contribution = SortContri;
VBPCAResult.FilteredData = FilteredData;
VBPCAResult.OriginTimeSeries = X;
VBPCAResult.SpatialResponsesMatrix = SpatialResponses'; % ATTENTION!!!
VBPCAResult.Components = Components; % Store sorted components
VBPCAResult.NoiseVariance = Result.V; % Store noise variance
end