function detection = offset_detect(ts, opts)
% Fits known MOM offsets first, then scans and classifies new candidates.

requiredFields = {'neu', 'mjd', 'decimalYear', 'station'};
[threadCleanup, threadsUsed, threadsAvailable] = ...
    limitComputationThreads(opts); %#ok<ASGLU>
for k = 1:numel(requiredFields)
    if ~isfield(ts, requiredFields{k})
        error('offset:MissingTimeSeriesField', ...
            'The time-series structure is missing field "%s".', requiredFields{k});
    end
end

Y = ts.neu;
[m, r] = size(Y);
if r ~= 3
    error('offset:InvalidComponents', 'Expected three NEU components.');
end
if numel(ts.mjd) ~= m || numel(ts.decimalYear) ~= m || ...
        any(~isfinite(Y), 'all') || any(~isfinite(ts.mjd(:))) || ...
        any(~isfinite(ts.decimalYear(:)))
    error('offset:InvalidTimeSeries', ...
        ['NEU, MJD, and decimal-year data must be finite and contain the ', ...
        'same number of observations.']);
end
if any(diff(ts.mjd(:)) <= 0)
    error('offset:UnorderedEpochs', 'MJD epochs must be strictly increasing.');
end

known = offset_knownEvents(ts);
knownIndices = known.indices(:).';
pointwiseThreshold = offset_chi2crit(opts.algorithm.alpha, r);
functionalA = offset_design(ts.decimalYear, opts.algorithm.harmonicCount, []);
functionalColumnCount = size(functionalA, 2);
baseA = offset_design(ts.decimalYear, ...
    opts.algorithm.harmonicCount, knownIndices);
[~, flickerEigenvalues, flickerEigenvectors] = ...
    offset_flicker(ts.mjd, false);
flickerCache = struct('eigenvalues', flickerEigenvalues, ...
    'eigenvectors', flickerEigenvectors);
baseNoise = offset_noise(Y, baseA, opts, ts.mjd, [], flickerCache);
baseFit = offset_fit(Y, baseA, baseNoise, opts);
baseNoiseSummary = compactNoise(baseNoise);

accepted = knownIndices;
acceptedPowers = nan(size(accepted));
acceptedThresholds = nan(size(accepted));
acceptedDetectionOrder = zeros(1, 0);
detectedIndices = zeros(1, 0);
diagnostics = repmat(emptyDiagnostic(), 0, 1);

A = baseA;
noise = baseNoise;
clear baseNoise
fit = baseFit;
iteration = 0;

while true
    iteration = iteration + 1;
    mask = offset_cand(m, accepted);
    [power, candidateIndices, offsetEstimates] = ...
        offset_power(Y, A, fit, noise, mask);

    fullPower = nan(m, 1);
    fullEstimates = nan(m, r);
    fullPower(candidateIndices) = power;
    fullEstimates(candidateIndices, :) = offsetEstimates;

    candidateIndex = NaN;
    candidateMjd = NaN;
    maxPower = NaN;
    candidateEstimate = nan(1, r);
    offsetMagnitudeMm = NaN;
    acceptedThisIteration = false;
    threshold = NaN;
    candidateCount = numel(candidateIndices);

    if fit.redundancy <= r
        terminationReason = 'insufficient-redundancy';
    elseif isempty(candidateIndices)
        terminationReason = 'no-identifiable-candidate';
    else
        threshold = offset_stationThreshold( ...
            opts.algorithm.alpha, r, candidateCount);
        [maxPower, localIndex] = max(power);
        candidateIndex = candidateIndices(localIndex);
        candidateMjd = ts.mjd(candidateIndex);
        candidateEstimate = offsetEstimates(localIndex, :);
        offsetMagnitudeMm = 1000 .* norm(candidateEstimate);
        acceptedThisIteration = isfinite(maxPower) && maxPower > threshold;
        if acceptedThisIteration
            terminationReason = 'accepted';
        else
            terminationReason = 'below-threshold';
        end
    end

    diagnostic = emptyDiagnostic();
    diagnostic.iteration = iteration;
    diagnostic.candidateIndex = candidateIndex;
    diagnostic.candidateMjd = candidateMjd;
    diagnostic.maxPower = maxPower;
    diagnostic.threshold = threshold;
    diagnostic.pointwiseThreshold = pointwiseThreshold;
    diagnostic.candidateCount = candidateCount;
    diagnostic.thresholdMethod = opts.algorithm.thresholdMethod;
    diagnostic.accepted = acceptedThisIteration;
    diagnostic.offsetEstimate = candidateEstimate;
    diagnostic.offsetMagnitudeMm = offsetMagnitudeMm;
    diagnostic.noiseModel = noise.model;
    diagnostic.whiteVariance = noise.whiteVariance;
    diagnostic.flickerVariance = noise.flickerVariance;
    diagnostic.componentCovariance = noise.componentCovariance;
    diagnostic.lsVceIterations = noise.lsVceIterations;
    diagnostic.lsVceConverged = noise.lsVceConverged;
    diagnostic.lsVceRelativeChange = noise.lsVceRelativeChange;
    diagnostic.lsVceNormalResidual = noise.lsVceNormalResidual;
    diagnostic.lsVceHistory = noise.lsVceHistory;
    diagnostic.warmStarted = noise.warmStarted;
    diagnostic.designRank = fit.designRank;
    diagnostic.redundancy = fit.redundancy;
    diagnostic.currentAcceptedIndices = accepted(:);
    diagnostic.candidateIndices = candidateIndices;
    diagnostic.candidateOffsetEstimates = fullEstimates;
    diagnostic.power = fullPower;
    diagnostic.terminationReason = terminationReason;
    diagnostics(end + 1, 1) = diagnostic; %#ok<AGROW>

    if ~acceptedThisIteration
        break;
    end

    accepted(end + 1) = candidateIndex; %#ok<AGROW>
    acceptedPowers(end + 1) = maxPower; %#ok<AGROW>
    acceptedThresholds(end + 1) = threshold; %#ok<AGROW>
    acceptedDetectionOrder(end + 1) = candidateIndex; %#ok<AGROW>
    detectedIndices(end + 1) = candidateIndex; %#ok<AGROW>
    [accepted, order] = sort(accepted);
    acceptedPowers = acceptedPowers(order);
    acceptedThresholds = acceptedThresholds(order);

    % Reuses the last variance ratio as the next initial value.
    warmStart = noise.componentVariances;
    clear noise
    A = offset_design(ts.decimalYear, opts.algorithm.harmonicCount, accepted);
    noise = offset_noise(Y, A, opts, ts.mjd, warmStart, flickerCache);
    fit = offset_fit(Y, A, noise, opts);
end

if isempty(accepted)
    offsetCoefficients = zeros(0, r);
else
    offsetCoefficients = fit.X(functionalColumnCount + 1:end, :);
end

knownEventMask = ismember(accepted(:), known.indices);
constraintMjd = nan(numel(accepted), 1);
[present, knownLocation] = ismember(accepted(:), known.indices);
constraintMjd(present) = known.requestedMjd(knownLocation(present));

detection = struct();
detection.station = ts.station;
detection.acceptedIndices = accepted(:);
detection.acceptedDetectionOrder = acceptedDetectionOrder(:);
detection.detectedIndices = sort(detectedIndices(:));
detection.knownEventIndices = known.indices(:);
detection.knownEventMjd = known.requestedMjd(:);
detection.knownEventAllMjd = known.allRequestedMjd(:);
detection.knownEventEffectiveMjd = known.effectiveMjd(:);
detection.knownEventMjdGroups = known.requestedMjdGroups;
detection.knownEventMask = knownEventMask;
detection.constraintMjd = constraintMjd;
detection.offsetCoefficients = offsetCoefficients;
detection.acceptedPowers = acceptedPowers(:);
detection.acceptedThresholds = acceptedThresholds(:);
detection.pointwiseThreshold = pointwiseThreshold;
detection.thresholdMethod = opts.algorithm.thresholdMethod;
detection.threshold = firstFiniteThreshold(diagnostics, pointwiseThreshold);
detection.alpha = opts.algorithm.alpha;
detection.threadsUsed = threadsUsed;
detection.threadsAvailable = threadsAvailable;
detection.baseFit = baseFit;
detection.baseNoise = baseNoiseSummary;
detection.detectionFit = fit;
detection.detectionNoise = compactNoise(noise);
detection.diagnostics = diagnostics;
detection.terminationReason = diagnostics(end).terminationReason;

classification = offset_classifyEvents(ts, detection, opts);
classificationFields = fieldnames(classification);
for k = 1:numel(classificationFields)
    fieldName = classificationFields{k};
    detection.(fieldName) = classification.(fieldName);
end

[detection.finalFit, finalNoise, detection.appliedIndices, ...
    detection.appliedCoefficients] = correctionModel( ...
    Y, ts, baseFit, functionalColumnCount, accepted, ...
    classification.applyCorrection, noise, baseNoiseSummary, opts, ...
    flickerCache);
detection.finalNoise = compactNoise(finalNoise);
detection.candidateCount = numel(detectedIndices);
detection.knownConstraintCount = numel(known.indices);
detection.knownConstraintSourceCount = numel(known.allRequestedMjd);
detection.collapsedKnownConstraintCount = known.collapsedCount;
detection.eventCount = numel(accepted);
detection.correctedOffsetCount = numel(detection.appliedIndices);
positivePairIds = unique(classification.pairId(classification.pairId > 0));
detection.transientPairCount = numel(positivePairIds);
detection.warnings = detectionWarnings(diagnostics, classification, known);
detection.offsetTable = offset_offsetTable(ts, detection);

end

function threshold = firstFiniteThreshold(diagnostics, fallback)
thresholds = [diagnostics.threshold];
first = find(isfinite(thresholds), 1, 'first');
if isempty(first)
    threshold = fallback;
else
    threshold = thresholds(first);
end
end

function [fit, noise, appliedIndices, appliedCoefficients] = ...
        correctionModel(Y, ts, baseFit, functionalColumnCount, ...
        accepted, applyCorrection, detectionNoise, baseNoiseSummary, ...
        opts, flickerCache)
appliedIndices = accepted(logical(applyCorrection(:).'));
if isempty(appliedIndices)
    fit = baseFit;
    noise = baseNoiseSummary;
    appliedCoefficients = zeros(0, size(Y, 2));
    return;
end

if isequal(appliedIndices, accepted)
    A = offset_design(ts.decimalYear, opts.algorithm.harmonicCount, accepted);
    fit = offset_fit(Y, A, detectionNoise, opts);
    noise = detectionNoise;
else
    A = offset_design(ts.decimalYear, ...
        opts.algorithm.harmonicCount, appliedIndices);
    warmStart = detectionNoise.componentVariances;
    noise = offset_noise(Y, A, opts, ts.mjd, warmStart, flickerCache);
    fit = offset_fit(Y, A, noise, opts);
end
appliedCoefficients = fit.X(functionalColumnCount + 1:end, :);
end

function [cleanup, threadsUsed, threadsAvailable] = limitComputationThreads(opts)
% Leaves CPU capacity for the desktop during dense calculations.
threadsUsed = NaN;
threadsAvailable = NaN;
cleanup = onCleanup(@() []);
if isfield(opts, 'runtime') && isstruct(opts.runtime) && ...
        isfield(opts.runtime, 'manageThreads') && ...
        ~logical(opts.runtime.manageThreads)
    try
        threadsAvailable = maxNumCompThreads;
        threadsUsed = 1;
    catch
        % Background workers manage their own thread capacity.
    end
    return;
end
previous = NaN;
try
    previous = maxNumCompThreads;
    threadsAvailable = previous;
    target = max(1, ceil(previous / 2));
    if target < previous
        maxNumCompThreads(target);
    end
    threadsUsed = maxNumCompThreads;
    cleanup = onCleanup(@() restoreComputationThreads(previous));
catch
    if isfinite(previous)
        restoreComputationThreads(previous);
    end
end
end

function restoreComputationThreads(threadCount)
try
    maxNumCompThreads(threadCount);
catch
    % Thread restoration errors do not replace the analysis result.
end
end

function summary = compactNoise(noise)
largeFields = {'QInv', 'precisionDesign', 'precisionResiduals', ...
    'stepQuadratic', 'fitCache'};
present = largeFields(isfield(noise, largeFields));
summary = rmfield(noise, present);
end

function warnings = detectionWarnings(diagnostics, classification, known)
warnings = cell(0, 1);
notConverged = find(~[diagnostics.lsVceConverged]);
if ~isempty(notConverged)
    warnings{end + 1, 1} = sprintf( ...
        'LS-VCE did not converge to tolerance at detection iteration(s): %s.', ...
        strjoin(string(notConverged), ', '));
end

if known.collapsedCount > 0
    warnings{end + 1, 1} = sprintf( ...
        ['%d known MOM offset date(s) were merged with another known ', ...
         'date because no retained common observation separates them; ', ...
         'each merged group is estimated as one net step.'], ...
        known.collapsedCount);
end

transientCount = nnz(classification.eventClass == "transient-pair");
if transientCount > 0
    warnings{end + 1, 1} = sprintf( ...
        ['%d detected epochs form %d short inverse pair(s) and were not ', ...
        'applied as permanent offset corrections.'], ...
        transientCount, transientCount ./ 2);
end
lowSupportCount = nnz(classification.eventClass == "low-support-candidate");
if lowSupportCount > 0
    warnings{end + 1, 1} = sprintf( ...
        ['%d detected epoch(s) lack the required observations or time span ', ...
        'on both sides and were not corrected automatically.'], ...
        lowSupportCount);
end

terminationReason = diagnostics(end).terminationReason;
switch terminationReason
    case 'insufficient-redundancy'
        warnings{end + 1, 1} = ...
            ['Detection stopped because adding another offset would leave ', ...
            'too few residual degrees of freedom to estimate the ', ...
            'multivariate component covariance.'];
    case 'no-identifiable-candidate'
        warnings{end + 1, 1} = ...
            ['Detection stopped because no remaining epoch supplied an ', ...
            'identifiable step column.'];
end
end

function diagnostic = emptyDiagnostic()
diagnostic = struct('iteration', 0, 'candidateIndex', NaN, ...
    'candidateMjd', NaN, 'maxPower', NaN, 'threshold', NaN, ...
    'pointwiseThreshold', NaN, 'candidateCount', 0, ...
    'thresholdMethod', '', ...
    'accepted', false, 'offsetEstimate', nan(1, 3), ...
    'offsetMagnitudeMm', NaN, 'noiseModel', '', ...
    'whiteVariance', NaN, 'flickerVariance', NaN, ...
    'componentCovariance', [], 'lsVceIterations', 0, ...
    'lsVceConverged', false, 'lsVceRelativeChange', NaN, ...
    'lsVceNormalResidual', NaN, 'lsVceHistory', [], ...
    'warmStarted', false, 'designRank', 0, 'redundancy', 0, ...
    'currentAcceptedIndices', [], 'candidateIndices', [], ...
    'candidateOffsetEstimates', [], 'power', [], ...
    'terminationReason', '');
end
