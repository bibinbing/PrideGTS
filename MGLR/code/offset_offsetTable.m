function offsetTable = offset_offsetTable(ts, detection)
% Creates an auditable table of all statistically accepted GLR candidates.

idx = detection.acceptedIndices(:);
n = numel(idx);
coeff = detection.offsetCoefficients;
if ~isequal(size(coeff), [n, 3])
    error('offset:OffsetTableAlignment', ...
        'Accepted epochs and three-component estimates are not aligned.');
end

dateValues = ts.time(idx);
dateValues.Format = 'yyyy-MM-dd';
dateString = string(dateValues);
magnitudeMm = 1000 .* sqrt(sum(coeff .^ 2, 2));
componentStdMm = candidateStandardErrors(detection, n);
power = alignedNumeric(detection, 'acceptedPowers', n, NaN);
threshold = alignedNumeric(detection, 'acceptedThresholds', n, ...
    detection.threshold);
pointwiseThreshold = repmat(detection.pointwiseThreshold, n, 1);

eventClass = alignedString(detection, 'eventClass', n, ...
    "persistent-candidate");
knownConstraint = alignedLogical(detection, 'knownEventMask', n, false);
constraintMjd = alignedNumeric(detection, 'constraintMjd', n, NaN);
qualityFlag = alignedString(detection, 'qualityFlag', n, "passed");
correctionReason = alignedString(detection, 'correctionReason', n, ...
    "accepted-persistent");
applyCorrection = alignedLogical(detection, 'applyCorrection', n, true);
pairId = alignedNumeric(detection, 'pairId', n, 0);
pairGapDays = alignedNumeric(detection, 'pairGapDays', n, NaN);
directionCosine = alignedNumeric(detection, 'directionCosine', n, NaN);
netToGrossRatio = alignedNumeric(detection, 'netToGrossRatio', n, NaN);
leftSupportObservations = alignedNumeric( ...
    detection, 'leftSupportObservations', n, NaN);
rightSupportObservations = alignedNumeric( ...
    detection, 'rightSupportObservations', n, NaN);
leftSupportDays = alignedNumeric(detection, 'leftSupportDays', n, NaN);
rightSupportDays = alignedNumeric(detection, 'rightSupportDays', n, NaN);

offsetTable = table(idx, ts.mjd(idx), ts.decimalYear(idx), dateString(:), ...
    coeff(:, 1), coeff(:, 2), coeff(:, 3), ...
    1000 .* coeff(:, 1), 1000 .* coeff(:, 2), 1000 .* coeff(:, 3), ...
    componentStdMm(:, 1), componentStdMm(:, 2), componentStdMm(:, 3), ...
    magnitudeMm, power, threshold, pointwiseThreshold, eventClass, ...
    knownConstraint, constraintMjd, qualityFlag, pairId, pairGapDays, ...
    directionCosine, netToGrossRatio, ...
    leftSupportObservations, rightSupportObservations, leftSupportDays, ...
    rightSupportDays, applyCorrection, correctionReason, ...
    'VariableNames', {'Index', 'MJD', 'DecimalYear', 'Date', ...
    'North_m', 'East_m', 'Up_m', 'North_mm', 'East_mm', 'Up_mm', ...
    'NorthStd_mm', 'EastStd_mm', 'UpStd_mm', 'Magnitude_mm', 'Power', ...
    'Threshold', 'PointwiseThreshold', 'EventClass', 'KnownConstraint', ...
    'ConstraintMJD', 'QualityFlag', 'PairID', 'PairGap_days', ...
    'DirectionCosine', 'NetToGrossRatio', ...
    'LeftSupportObservations', 'RightSupportObservations', ...
    'LeftSupport_days', 'RightSupport_days', 'ApplyCorrection', ...
    'CorrectionReason'});

end

function componentStdMm = candidateStandardErrors(detection, n)
componentStdMm = nan(n, 3);
if n == 0
    return;
end
if isfield(detection, 'detectionFit')
    fit = detection.detectionFit;
elseif isfield(detection, 'finalFit')
    fit = detection.finalFit;
else
    return;
end
if ~isfield(fit, 'Sigma')
    return;
end

normalCovariance = [];
if isfield(fit, 'normalInverse')
    normalCovariance = fit.normalInverse;
elseif isfield(fit, 'normalSolver')
    columnCount = size(fit.normalSolver.R, 1);
    normalCovariance = offset_normalSolve([], eye(columnCount), [], ...
        fit.normalSolver);
end
if isempty(normalCovariance) || size(normalCovariance, 1) < n
    return;
end
coefficientColumns = size(normalCovariance, 1) - n + (1:n);
temporalVariance = diag(normalCovariance(coefficientColumns, coefficientColumns));
componentVariance = max(real(diag(fit.Sigma)), 0).';
componentStdMm = 1000 .* sqrt(max(real(temporalVariance), 0) .* componentVariance);
end

function value = alignedNumeric(source, fieldName, n, fallback)
value = repmat(fallback, n, 1);
if isfield(source, fieldName)
    proposed = source.(fieldName);
    if isnumeric(proposed) && numel(proposed) == n
        value = proposed(:);
    end
end
end

function value = alignedString(source, fieldName, n, fallback)
value = repmat(string(fallback), n, 1);
if isfield(source, fieldName)
    proposed = string(source.(fieldName));
    if numel(proposed) == n
        value = proposed(:);
    end
end
end

function value = alignedLogical(source, fieldName, n, fallback)
value = repmat(logical(fallback), n, 1);
if isfield(source, fieldName)
    proposed = source.(fieldName);
    if (islogical(proposed) || isnumeric(proposed)) && numel(proposed) == n
        value = logical(proposed(:));
    end
end
end
