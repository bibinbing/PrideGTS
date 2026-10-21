function fit = offset_fit(Y, A, noise, opts)
% Fits the multivariate functional model by generalized least squares.

[m, r] = size(Y);
if size(A, 1) ~= m || isempty(A) || any(~isfinite(A), 'all') || ...
        any(~isfinite(Y), 'all')
    error('offset:InvalidFitInput', ...
        'The observations and design matrix must be finite with equal row counts.');
end

p = size(A, 2);
designRank = rank(A);
redundancy = m - designRank;
if designRank ~= p || redundancy < r
    error('offset:UnidentifiableFunctionalModel', ...
        ['The functional model must have full column rank and at least one ', ...
        'residual degree of freedom per component to estimate the %d-by-%d ', ...
        'component covariance (rank %d, columns %d, observations %d).'], ...
        r, r, designRank, p, m);
end

useCache = isfield(noise, 'fitCache') && isstruct(noise.fitCache) && ...
    isfield(noise.fitCache, 'design') && isequal(noise.fitCache.design, A) && ...
    isfield(noise.fitCache, 'observations') && ...
    isequal(noise.fitCache.observations, Y) && ...
    isfield(noise.fitCache, 'X') && isequal(size(noise.fitCache.X), [p, r]) && ...
    isfield(noise.fitCache, 'residuals') && ...
    isequal(size(noise.fitCache.residuals), [m, r]);
if useCache
    cache = noise.fitCache;
    X = cache.X;
    residuals = cache.residuals;
    Sigma = cache.componentCovariance;
    normal = cache.normal;
    normalSolver = cache.normalSolver;
else
    QInv = noise.QInv;
    if ~isequal(size(QInv), [m, m]) || any(~isfinite(QInv), 'all')
        error('offset:InvalidTemporalPrecision', ...
            ['Noise information must contain either a matching exact fit ', ...
            'cache or a finite m-by-m temporal precision matrix.']);
    end

    QInvA = QInv * A;
    normal = A.' * QInvA;
    rhs = A.' * (QInv * Y);
    [X, normalSolver] = offset_normalSolve(normal, rhs, ...
        opts.algorithm.maxConditionNumber);
    residuals = Y - A * X;
    Sigma = (residuals.' * (QInv * residuals)) ./ redundancy;
    Sigma = regularizeComponentCovariance(Sigma, opts);
end

fit = struct();
fit.X = X;
fit.residuals = residuals;
fit.Sigma = Sigma;
fit.redundancy = redundancy;
fit.designRank = designRank;
fit.design = A;
fit.normal = normal;
fit.normalSolver = normalSolver;
fit.normalRcond = normalSolver.reciprocalCondition;
fit.maxConditionNumber = opts.algorithm.maxConditionNumber;
fit.componentCount = r;

end

function Sigma = regularizeComponentCovariance(Sigma, opts)
Sigma = (Sigma + Sigma.') ./ 2;
scale = max(abs(diag(Sigma)));
if isempty(scale) || ~isfinite(scale) || scale <= 0
    error('offset:InvalidComponentCovariance', ...
        'The residual component covariance is not positive and finite.');
end
jitter = max(eps(scale), opts.noise.ridge .* scale);
[~, flag] = chol(Sigma);
while flag ~= 0
    Sigma = Sigma + eye(size(Sigma)) .* jitter;
    jitter = jitter .* 10;
    [~, flag] = chol(Sigma);
end

end
