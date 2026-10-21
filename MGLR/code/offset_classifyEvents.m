function classification = offset_classifyEvents(ts, detection, opts)
% Classifies accepted candidates before coordinate correction.

indices = detection.acceptedIndices(:);
coefficients = detection.offsetCoefficients;
n = numel(indices);
if ~isequal(size(coefficients), [n, 3])
    error('offset:ClassificationAlignment', ...
        'Accepted epochs and three-component offset estimates are not aligned.');
end

classification = emptyClassification(n);
if n == 0
    return;
end

classification.eventClass(:) = "persistent-candidate";
classification.qualityFlag(:) = "passed";
classification.applyCorrection(:) = true;
classification.correctionReason(:) = "accepted-persistent";
knownMask = alignedKnownMask(detection, indices);

[leftCount, rightCount, leftDays, rightDays] = ...
    segmentSupport(ts.mjd(:), indices);
classification.leftSupportObservations = leftCount;
classification.rightSupportObservations = rightCount;
classification.leftSupportDays = leftDays;
classification.rightSupportDays = rightDays;

lowSupport = ~knownMask & (leftCount < opts.quality.minSideObservations | ...
    rightCount < opts.quality.minSideObservations | ...
    leftDays < opts.quality.minSideSpanDays | ...
    rightDays < opts.quality.minSideSpanDays);
classification.eventClass(lowSupport) = "low-support-candidate";
classification.qualityFlag(lowSupport) = "low-support";
classification.applyCorrection(lowSupport) = false;
classification.correctionReason(lowSupport) = "not-applied-low-support";

[pairCandidates, pairGap, directionCosine, netToGross] = ...
    inversePairMetrics(ts.mjd(:), indices, coefficients, opts);
knownAdjacent = knownMask(1:end - 1) | knownMask(2:end);
pairCandidates(knownAdjacent) = false;
pairGap(knownAdjacent) = NaN;
directionCosine(knownAdjacent) = NaN;
netToGross(knownAdjacent) = NaN;
classification.pairGapDays(1:end - 1) = pairGap;
classification.directionCosine(1:end - 1) = directionCosine;
classification.netToGrossRatio(1:end - 1) = netToGross;

% Selects non-overlapping pairs by cancellation strength.
candidateStarts = find(pairCandidates);
if ~isempty(candidateStarts)
    [~, localOrder] = sort(netToGross(candidateStarts), 'ascend');
    candidateStarts = candidateStarts(localOrder);
end
used = false(n, 1);
pairId = 0;
for k = 1:numel(candidateStarts)
    first = candidateStarts(k);
    second = first + 1;
    if used(first) || used(second)
        continue;
    end
    pairId = pairId + 1;
    used([first, second]) = true;
    classification.pairId([first, second]) = pairId;
    classification.pairGapDays([first, second]) = pairGap(first);
    classification.directionCosine([first, second]) = directionCosine(first);
    classification.netToGrossRatio([first, second]) = netToGross(first);
    classification.eventClass([first, second]) = "transient-pair";
    classification.qualityFlag([first, second]) = ...
        appendFlag(classification.qualityFlag([first, second]), ...
        "short-inverse-pair");
    classification.applyCorrection([first, second]) = false;
    classification.correctionReason([first, second]) = ...
        "not-applied-transient-pair";
end

classification.eventClass(knownMask) = "known-constrained";
classification.qualityFlag(knownMask) = "known-event";
classification.applyCorrection(knownMask) = true;
classification.correctionReason(knownMask) = "applied-known-constraint";

end

function knownMask = alignedKnownMask(detection, indices)
knownMask = false(numel(indices), 1);
if isfield(detection, 'knownEventMask')
    proposed = detection.knownEventMask;
    if numel(proposed) ~= numel(indices)
        error('offset:ClassificationAlignment', ...
            'Known-event flags and accepted epochs are not aligned.');
    end
    knownMask = logical(proposed(:));
elseif isfield(detection, 'knownEventIndices')
    knownMask = ismember(indices, detection.knownEventIndices(:));
end
end

function classification = emptyClassification(n)
classification = struct();
classification.eventClass = strings(n, 1);
classification.qualityFlag = strings(n, 1);
classification.applyCorrection = false(n, 1);
classification.correctionReason = strings(n, 1);
classification.pairId = zeros(n, 1);
classification.pairGapDays = nan(n, 1);
classification.directionCosine = nan(n, 1);
classification.netToGrossRatio = nan(n, 1);
classification.leftSupportObservations = zeros(n, 1);
classification.rightSupportObservations = zeros(n, 1);
classification.leftSupportDays = zeros(n, 1);
classification.rightSupportDays = zeros(n, 1);
end

function [leftCount, rightCount, leftDays, rightDays] = ...
        segmentSupport(mjd, indices)
n = numel(indices);
leftCount = zeros(n, 1);
rightCount = zeros(n, 1);
leftDays = zeros(n, 1);
rightDays = zeros(n, 1);
for i = 1:n
    if i == 1
        leftBoundary = 1;
    else
        leftBoundary = indices(i - 1);
    end
    if i == n
        rightBoundary = numel(mjd);
    else
        rightBoundary = indices(i + 1) - 1;
    end
    leftCount(i) = indices(i) - leftBoundary;
    rightCount(i) = rightBoundary - indices(i) + 1;
    leftDays(i) = mjd(indices(i)) - mjd(leftBoundary);
    rightDays(i) = mjd(rightBoundary) - mjd(indices(i));
end
end

function [qualifies, gapDays, directionCosine, netToGross] = ...
        inversePairMetrics(mjd, indices, coefficients, opts)
n = numel(indices);
if n < 2
    qualifies = false(0, 1);
    gapDays = zeros(0, 1);
    directionCosine = zeros(0, 1);
    netToGross = zeros(0, 1);
    return;
end

first = coefficients(1:end - 1, :);
second = coefficients(2:end, :);
firstNorm = sqrt(sum(first .^ 2, 2));
secondNorm = sqrt(sum(second .^ 2, 2));
gross = firstNorm + secondNorm;
gapDays = mjd(indices(2:end)) - mjd(indices(1:end - 1));
directionCosine = sum(first .* second, 2) ./ max(firstNorm .* secondNorm, realmin);
netToGross = sqrt(sum((first + second) .^ 2, 2)) ./ max(gross, realmin);
finite = isfinite(gapDays) & isfinite(directionCosine) & ...
    isfinite(netToGross) & firstNorm > 0 & secondNorm > 0;
qualifies = finite & gapDays <= opts.quality.transientMaxGapDays & ...
    directionCosine <= opts.quality.transientMaxDirectionCosine & ...
    netToGross <= opts.quality.transientMaxNetToGrossRatio;
end

function flags = appendFlag(flags, newFlag)
flags = string(flags);
for i = 1:numel(flags)
    if flags(i) == "" || flags(i) == "passed"
        flags(i) = newFlag;
    elseif ~contains(flags(i), newFlag)
        flags(i) = flags(i) + ";" + newFlag;
    end
end
end
