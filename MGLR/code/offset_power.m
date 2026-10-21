function [power, candidateIndices, offsetEstimates] = offset_power(Y, A, fit, noise, mask)
% Evaluates step candidates without forming a dense precision matrix.

m = size(Y, 1);
r = size(Y, 2);
if numel(mask) ~= m
    error('offset:CandidateMaskSize', ...
        'The candidate mask must contain one value per observation.');
end

requestedIndices = find(mask(:));
if isempty(requestedIndices) || fit.redundancy <= r
    power = zeros(0, 1);
    candidateIndices = zeros(0, 1);
    offsetEstimates = zeros(0, r);
    return;
end

[precisionDesign, precisionResiduals, direct] = ...
    precisionSufficientStatistics(A, fit, noise);
stepDesign = flipud(cumsum(flipud(precisionDesign), 1));
stepResiduals = flipud(cumsum(flipud(precisionResiduals), 1));

v = stepDesign(requestedIndices, :).';
projectedV = offset_normalSolve([], v, [], fit.normalSolver);
projection = sum(v .* projectedV, 1).';
direct = direct(requestedIndices);
denominator = direct - projection;

% Keeps candidates that add an independent step column.
scale = max(abs(direct), abs(projection));
scale = max(scale, realmin);
denominatorTolerance = max(100 .* eps(scale), ...
    scale ./ fit.maxConditionNumber);
valid = isfinite(denominator) & direct > 0 & ...
    denominator > denominatorTolerance;

candidateIndices = requestedIndices(valid);
denominator = denominator(valid);
numerator = stepResiduals(candidateIndices, :).';
sigmaWeightedNumerator = fit.Sigma \ numerator;
power = sum(numerator .* sigmaWeightedNumerator, 1).' ./ denominator;
power = max(0, real(power));
offsetEstimates = (numerator ./ denominator.').';

finite = isfinite(power) & all(isfinite(offsetEstimates), 2);
candidateIndices = candidateIndices(finite);
power = power(finite);
offsetEstimates = offsetEstimates(finite, :);

end

function [precisionDesign, precisionResiduals, direct] = ...
        precisionSufficientStatistics(A, fit, noise)
m = size(A, 1);
if isfield(noise, 'precisionDesign') && ...
        isequal(size(noise.precisionDesign), size(A)) && ...
        isfield(noise, 'precisionResiduals') && ...
        isequal(size(noise.precisionResiduals), size(fit.residuals)) && ...
        isfield(noise, 'stepQuadratic') && numel(noise.stepQuadratic) == m
    precisionDesign = noise.precisionDesign;
    precisionResiduals = noise.precisionResiduals;
    direct = noise.stepQuadratic(:);
    return;
end

% Supports a legacy dense cache when supplied by a caller.
if ~isfield(noise, 'QInv') || ~isequal(size(noise.QInv), [m, m])
    error('offset:MissingPrecisionStatistics', ...
        'Noise information lacks the exact precision sufficient statistics.');
end
QInv = noise.QInv;
precisionDesign = QInv * A;
precisionResiduals = QInv * fit.residuals;
columnSums = cumsum(QInv(:, end:-1:1), 2);
columnSums = columnSums(:, end:-1:1);
direct = zeros(m, 1);
for j = 1:m
    direct(j) = sum(columnSums(j:end, j));
end
end
