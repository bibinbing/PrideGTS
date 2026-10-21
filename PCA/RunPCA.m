function PCAResult = RunPCA(X, NumReconComp)
% RunPCA - PCA-based spatiotemporal filtering without extra toolboxes
%
% Inputs:
%   X            - m × n residual time-series matrix
%                  m: epochs/days, n: stations
%   NumReconComp - Number of principal components used to reconstruct CME
%
% Output:
%   PCAResult    - Structure containing the PCA and filtering results

% Input validation
if ~isnumeric(X) || ~ismatrix(X) || ~isreal(X) || isempty(X)
    error('X must be a nonempty real numeric matrix.');
end

if any(~isfinite(X(:)))
    error('X contains NaN or Inf. Please preprocess missing data first.');
end

X = double(X);
[numEpochs, numStations] = size(X);

if numEpochs < 2
    error('X must contain at least two observations.');
end

numComponents = min(numEpochs - 1, numStations);

if ~isnumeric(NumReconComp) || ~isscalar(NumReconComp) || ...
        ~isfinite(NumReconComp) || NumReconComp < 0 || ...
        NumReconComp ~= floor(NumReconComp) || ...
        NumReconComp > numComponents
    error('NumReconComp must be an integer between 0 and %d.', ...
          numComponents);
end

% Center the observations
dataMean = mean(X, 1);
XCentered = bsxfun(@minus, X, dataMean);

% PCA using singular value decomposition
% XCentered = U*S*V'
% Columns of V are the principal-component loading vectors.
[U, S, V] = svd(XCentered, 'econ');

coeff = V(:, 1:numComponents);
score = U(:, 1:numComponents) * ...
        S(1:numComponents, 1:numComponents);

singularValues = diag(S);
singularValues = singularValues(1:numComponents);

% Eigenvalues of the sample covariance matrix
latent = singularValues.^2 / (numEpochs - 1);

% Percentage variance explained by every component
totalVariance = sum(latent);

if totalVariance > 0
    explained = 100 * latent / totalVariance;
else
    explained = zeros(size(latent));
end

% Reconstruct the common-mode error
if NumReconComp == 0
    cme = zeros(size(X));
else
    cme = score(:, 1:NumReconComp) * ...
          coeff(:, 1:NumReconComp)';
end

% Remove the reconstructed common-mode error
% The CME is reconstructed from centered data, so the station means must
% only be restored once.
XFilteredCentered = XCentered - cme;
XFiltered = bsxfun(@plus, XFilteredCentered, dataMean);

% Store results
PCAResult.OriginTimeSeries       = X;
PCAResult.FilteredData           = XFiltered;
PCAResult.CME                    = cme;
PCAResult.Mean                   = dataMean;
PCAResult.Contribution           = explained;
PCAResult.Eigenvalues            = latent;
PCAResult.SpatialResponsesMatrix = coeff;
PCAResult.Components             = score;

end