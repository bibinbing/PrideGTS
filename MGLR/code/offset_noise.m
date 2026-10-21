function noise = offset_noise(Y, A, opts, epochDays, warmStart, flickerCache)
% Reuses one eigendecomposition while updating the noise weights.

[m, r] = size(Y);
if nargin < 4 || isempty(epochDays)
    epochDays = (0:m - 1).';
end
if nargin < 5
    warmStart = [];
end
if nargin < 6
    flickerCache = [];
end
if size(A, 1) ~= m || any(~isfinite(Y), 'all') || any(~isfinite(A), 'all')
    error('offset:InvalidNoiseInput', ...
        'The observations and design matrix must be finite with equal row counts.');
end
if m - size(A, 2) < r || rank(A) ~= size(A, 2)
    error('offset:UnidentifiableNoiseModel', ...
        ['LS-VCE requires a full-rank design with at least %d residual ', ...
        'degrees of freedom to estimate the component covariance.'], r);
end

model = normalizeNoiseModel(opts.noise.model);
iterationLimit = validatePositiveInteger(opts.noise.lsVceIterations, ...
    'opts.noise.lsVceIterations');
tolerance = opts.noise.lsVceTolerance;
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error('offset:InvalidLsVceTolerance', ...
        'opts.noise.lsVceTolerance must be a finite positive scalar.');
end

if isempty(flickerCache)
    [~, flickerEigenvalues, flickerEigenvectors] = offset_flicker(epochDays);
else
    if ~isstruct(flickerCache) || ...
            ~isfield(flickerCache, 'eigenvalues') || ...
            ~isfield(flickerCache, 'eigenvectors')
        error('offset:InvalidFlickerCache', ...
            'The flicker cache must contain eigenvalues and eigenvectors.');
    end
    flickerEigenvalues = flickerCache.eigenvalues;
    flickerEigenvectors = flickerCache.eigenvectors;
    if numel(flickerEigenvalues) ~= m || ...
            ~isequal(size(flickerEigenvectors), [m, m]) || ...
            any(~isfinite(flickerEigenvalues)) || ...
            any(~isfinite(flickerEigenvectors), 'all')
        error('offset:InvalidFlickerCache', ...
            'The flicker cache does not match the observation count.');
    end
    flickerEigenvalues = flickerEigenvalues(:);
end
meanFlickerDiagonal = mean(flickerEigenvalues);
[theta, warmStarted] = initialTemporalComponents( ...
    warmStart, meanFlickerDiagonal, opts);
history = repmat(emptyHistoryEntry(), iterationLimit, 1);
converged = false;

for iter = 1:iterationLimit
    thetaStart = theta;
    [factor, temporalJitter, temporalRcond] = ...
        temporalFactor(flickerEigenvectors, flickerEigenvalues, theta, opts);
    state = glsStateFromFactor(Y, A, factor, opts);
    effectiveWhite = theta(1) + temporalJitter;
    [thetaEstimate, varianceNormal, varianceRhs, normalResidual, traceMethod] = ...
        estimateTemporalComponents(state, factor, A, ...
        effectiveWhite, theta(2), r, opts);
    theta = normalizeTemporalComponents( ...
        thetaEstimate, meanFlickerDiagonal, opts);
    relativeChange = norm(theta - thetaStart, inf) / ...
        max(norm(thetaStart, inf), opts.noise.minVariance);

    history(iter).iteration = iter;
    history(iter).thetaStart = thetaStart;
    history(iter).thetaEstimate = thetaEstimate;
    history(iter).thetaNormalized = theta;
    history(iter).relativeChange = relativeChange;
    history(iter).normalEquationResidual = normalResidual;
    history(iter).varianceNormalMatrix = varianceNormal;
    history(iter).varianceRightHandSide = varianceRhs;
    history(iter).componentCovariance = state.componentCovariance;
    history(iter).weightedResidualSumSquares = ...
        sum(state.residuals .* state.weightedResiduals, 1).';
    history(iter).temporalJitter = temporalJitter;
    history(iter).temporalReciprocalCondition = temporalRcond;
    history(iter).traceMethod = traceMethod;

    if relativeChange <= tolerance
        converged = true;
        break;
    end
end

iterationsUsed = iter;
history = history(1:iterationsUsed);

% Builds the final fit cache used by the candidate scan.
[factor, temporalJitter, temporalRcond] = ...
    temporalFactor(flickerEigenvectors, flickerEigenvalues, theta, opts);
finalState = glsStateFromFactor(Y, A, factor, opts);
stepQuadratic = stepQuadraticFromFactor(factor);

noise = struct();
noise.isDiagonal = false;
noise.weights = [];
noise.QInv = [];
noise.precisionDesign = finalState.QInvA;
noise.precisionResiduals = finalState.weightedResiduals;
noise.stepQuadratic = stepQuadratic;
noise.precisionRepresentation = 'factorized-exact';
noise.whiteVariance = theta(1);
noise.flickerVariance = theta(2);
noise.model = [model '-lsvce'];
noise.componentVariances = theta(:);
noise.componentCovariance = finalState.componentCovariance;
noise.lsVceIterations = iterationsUsed;
noise.lsVceConverged = converged;
noise.lsVceRelativeChange = history(end).relativeChange;
noise.lsVceNormalResidual = history(end).normalEquationResidual;
noise.lsVceTolerance = tolerance;
noise.lsVceHistory = history;
noise.lsVceTermination = ternary(converged, 'tolerance', 'iteration-limit');
noise.warmStarted = warmStarted;
noise.temporalJitter = temporalJitter;
noise.temporalReciprocalCondition = temporalRcond;
noise.fitCache = finalState;
noise.fitCache.design = A;
noise.fitCache.observations = Y;

end

function model = normalizeNoiseModel(model)
model = char(lower(strtrim(string(model))));
switch model
    case {'white-flicker', 'flicker-white'}
        model = 'white-flicker';
    otherwise
        error('offset:UnknownNoiseModel', ...
            ['The detector requires the multivariate white-plus-flicker ', ...
            'LS-VCE model; received "%s".'], model);
end
end

function value = validatePositiveInteger(value, label)
if ~isscalar(value) || ~isfinite(value) || value < 1 || value ~= fix(value)
    error('offset:InvalidPositiveInteger', '%s must be a positive integer.', label);
end
end

function [theta, warmStarted] = initialTemporalComponents( ...
        warmStart, meanFlickerDiagonal, opts)
theta = [1; 1];
warmStarted = false;
if ~isempty(warmStart)
    if isstruct(warmStart) && isfield(warmStart, 'componentVariances')
        proposed = warmStart.componentVariances;
    elseif isnumeric(warmStart)
        proposed = warmStart;
    else
        proposed = [];
    end
    if numel(proposed) == 2 && all(isfinite(proposed(:))) && all(proposed(:) > 0)
        theta = proposed(:);
        warmStarted = true;
    end
end
theta = normalizeTemporalComponents(theta, meanFlickerDiagonal, opts);
end

function theta = normalizeTemporalComponents(theta, meanFlickerDiagonal, opts)
theta = max(real(theta(:)), opts.noise.minVariance);
scale = theta(1) + theta(2) * meanFlickerDiagonal;
if ~isfinite(scale) || scale <= 0
    error('offset:InvalidTemporalVariance', ...
        'LS-VCE produced an invalid temporal covariance scale.');
end
theta = theta / scale;
theta = max(theta, opts.noise.minVariance);
scale = theta(1) + theta(2) * meanFlickerDiagonal;
theta = theta / scale;
end

function [factor, jitter, reciprocalCondition] = ...
        temporalFactor(flickerEigenvectors, flickerEigenvalues, theta, opts)
covarianceEigenvalues = theta(1) + theta(2) .* flickerEigenvalues;
jitter = 0;
if any(~isfinite(covarianceEigenvalues))
    error('offset:InvalidTemporalSpectrum', ...
        'The temporal covariance spectrum must be finite.');
end
if any(covarianceEigenvalues <= 0)
    scale = max(mean(abs(covarianceEigenvalues)), 1);
    trialJitter = max(eps(scale), opts.noise.ridge * scale);
    while any(covarianceEigenvalues <= 0)
        covarianceEigenvalues = covarianceEigenvalues + trialJitter;
        jitter = jitter + trialJitter;
        trialJitter = trialJitter * 10;
        if jitter > scale
            error('offset:TemporalCovarianceNotPositiveDefinite', ...
                'The LS-VCE temporal covariance could not be factorized.');
        end
    end
end
precisionEigenvalues = 1 ./ covarianceEigenvalues;
reciprocalCondition = min(covarianceEigenvalues) / max(covarianceEigenvalues);
factor = struct( ...
    'eigenvectors', flickerEigenvectors, ...
    'flickerEigenvalues', flickerEigenvalues, ...
    'precisionEigenvalues', precisionEigenvalues, ...
    'sqrtPrecisionEigenvalues', sqrt(precisionEigenvalues));
end

function X = solveTemporal(factor, B)
spectralB = factor.eigenvectors.' * B;
X = factor.eigenvectors * (factor.precisionEigenvalues .* spectralB);
end

function FB = applyFlicker(factor, B)
spectralB = factor.eigenvectors.' * B;
FB = factor.eigenvectors * (factor.flickerEigenvalues .* spectralB);
end

function stepQuadratic = stepQuadraticFromFactor(factor)
% Computes each step term without storing a precision matrix.
m = size(factor.eigenvectors, 1);
stepQuadratic = zeros(m, 1);
spectralStep = zeros(1, m);
precision = factor.precisionEigenvalues.';
for j = m:-1:1
    spectralStep = spectralStep + factor.eigenvectors(j, :);
    stepQuadratic(j) = sum((spectralStep .^ 2) .* precision);
end
end

function state = glsStateFromFactor(Y, A, factor, opts)
columnCount = size(A, 2);
spectralData = factor.eigenvectors.' * [A, Y];
whitenedData = factor.sqrtPrecisionEigenvalues .* spectralData;
whitenedA = whitenedData(:, 1:columnCount);
whitenedY = whitenedData(:, columnCount + 1:end);
normal = whitenedA.' * whitenedA;
rhs = whitenedA.' * whitenedY;
[X, normalSolver] = offset_normalSolve(normal, rhs, ...
    opts.algorithm.maxConditionNumber);
residuals = Y - A * X;
whitenedResiduals = whitenedY - whitenedA * X;
spectralResiduals = spectralData(:, columnCount + 1:end) - ...
    spectralData(:, 1:columnCount) * X;
precisionData = factor.eigenvectors * ...
    (factor.precisionEigenvalues .* ...
    [spectralData(:, 1:columnCount), spectralResiduals]);
QInvA = precisionData(:, 1:columnCount);
weightedResiduals = precisionData(:, columnCount + 1:end);
redundancy = size(Y, 1) - size(A, 2);
componentCovariance = (whitenedResiduals.' * whitenedResiduals) ./ redundancy;
componentCovariance = regularizeComponentCovariance(componentCovariance, opts);

state = struct('X', X, 'residuals', residuals, ...
    'weightedResiduals', weightedResiduals, ...
    'componentCovariance', componentCovariance, 'QInvA', QInvA, ...
    'normal', normal, 'normalSolver', normalSolver);
end

function Sigma = regularizeComponentCovariance(Sigma, opts)
Sigma = (Sigma + Sigma.') ./ 2;
scale = max(abs(diag(Sigma)));
if isempty(scale) || ~isfinite(scale) || scale <= 0
    error('offset:InvalidComponentCovariance', ...
        'The residual component covariance is not positive and finite.');
end
jitter = max(eps(scale), opts.noise.ridge * scale);
[~, flag] = chol(Sigma);
while flag ~= 0
    Sigma = Sigma + eye(size(Sigma)) .* jitter;
    jitter = jitter * 10;
    [~, flag] = chol(Sigma);
end
end

function [theta, G, ell, normalResidual, traceMethod] = ...
        estimateTemporalComponents(state, factor, A, ...
        whiteLevel, flickerLevel, r, opts)
% Evaluates the variance-component trace terms in the shared basis.

B = state.QInvA;
flickerEigenvalues = factor.flickerEigenvalues;
precisionEigenvalues = factor.precisionEigenvalues;
weightedSquare = precisionEigenvalues .^ 2;
baseTrace = [sum(weightedSquare), ...
    sum(flickerEigenvalues .* weightedSquare); ...
    sum(flickerEigenvalues .* weightedSquare), ...
    sum((flickerEigenvalues .^ 2) .* weightedSquare)];

KB = solveTemporal(factor, B);
useIdentity = flickerLevel > sqrt(eps) * max(whiteLevel, 1);
if useIdentity
    FB = (A - whiteLevel .* B) ./ flickerLevel;
    KFB = (B - whiteLevel .* KB) ./ flickerLevel;
    FW = (state.residuals - whiteLevel .* state.weightedResiduals) ./ ...
        flickerLevel;
    traceMethod = 'spectral-trace-identity';
else
    FB = applyFlicker(factor, B);
    KFB = solveTemporal(factor, FB);
    FW = applyFlicker(factor, state.weightedResiduals);
    traceMethod = 'spectral-trace-direct-boundary';
end

D = cell(2, 1);
D{1} = B.' * B;
D{2} = B.' * FB;
H = cell(2, 2);
H{1, 1} = B.' * KB;
H{1, 2} = B.' * KFB;
H{2, 1} = H{1, 2}.';
H{2, 2} = FB.' * KFB;

G = zeros(2, 2);
CD = cell(2, 1);
for i = 1:2
    CD{i} = offset_normalSolve([], D{i}, [], state.normalSolver);
end
for i = 1:2
    for j = i:2
        CH = offset_normalSolve([], H{i, j}, [], state.normalSolver);
        projectedTrace = baseTrace(i, j) - 2 * trace(CH) + ...
            trace(CD{i} * CD{j});
        G(i, j) = r * projectedTrace;
        G(j, i) = G(i, j);
    end
end
G = (G + G.') ./ 2;

sigmaWeightedW = (state.componentCovariance \ state.weightedResiduals.').';
ell = [sum(sigmaWeightedW .* state.weightedResiduals, 'all'); ...
    sum(sigmaWeightedW .* FW, 'all')];

[theta, normalResidual] = constrainedVarianceSolve(G, ell, ...
    opts.noise.minVariance);
end

function [theta, projectedResidual] = constrainedVarianceSolve(G, ell, lowerBound)
% Solves the two constrained variance components.
G = (G + G.') ./ 2;
candidates = zeros(2, 4);
candidateCount = 0;

if rcond(G) > eps
    interior = G \ ell;
    if all(isfinite(interior)) && all(interior >= lowerBound)
        candidateCount = candidateCount + 1;
        candidates(:, candidateCount) = interior;
    end
end

theta2 = (ell(2) - G(2, 1) * lowerBound) / max(G(2, 2), realmin);
candidateCount = candidateCount + 1;
candidates(:, candidateCount) = [lowerBound; max(theta2, lowerBound)];
theta1 = (ell(1) - G(1, 2) * lowerBound) / max(G(1, 1), realmin);
candidateCount = candidateCount + 1;
candidates(:, candidateCount) = [max(theta1, lowerBound); lowerBound];
candidateCount = candidateCount + 1;
candidates(:, candidateCount) = [lowerBound; lowerBound];
candidates = candidates(:, 1:candidateCount);

objective = inf(1, size(candidates, 2));
for k = 1:size(candidates, 2)
    candidate = candidates(:, k);
    if all(isfinite(candidate))
        objective(k) = 0.5 * candidate.' * G * candidate - ell.' * candidate;
    end
end
[~, best] = min(objective);
theta = candidates(:, best);

gradient = G * theta - ell;
active = theta <= lowerBound * (1 + 10 * eps);
gradient(active) = min(gradient(active), 0);
projectedResidual = norm(gradient) / max(norm(ell), eps);
end

function entry = emptyHistoryEntry()
entry = struct('iteration', 0, 'thetaStart', [], 'thetaEstimate', [], ...
    'thetaNormalized', [], 'relativeChange', NaN, ...
    'normalEquationResidual', NaN, 'varianceNormalMatrix', [], ...
    'varianceRightHandSide', [], 'componentCovariance', [], ...
    'weightedResidualSumSquares', [], 'temporalJitter', NaN, ...
    'temporalReciprocalCondition', NaN, 'traceMethod', '');
end

function out = ternary(condition, trueValue, falseValue)
if condition
    out = trueValue;
else
    out = falseValue;
end
end
