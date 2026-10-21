function [ts, quality] = offset_cleanSeries(ts)
% Removes invalid MOM epochs and keeps the three components aligned.

requiredFields = {'mjd', 'time', 'decimalYear', 'neu'};
for k = 1:numel(requiredFields)
    if ~isfield(ts, requiredFields{k})
        error('offset:MissingTimeSeriesField', ...
            'The time series is missing the required field "%s".', ...
            requiredFields{k});
    end
end

n = numel(ts.mjd);
if numel(ts.time) ~= n || numel(ts.decimalYear) ~= n || ...
        size(ts.neu, 1) ~= n || size(ts.neu, 2) ~= 3
    error('offset:RowAlignment', ...
        'MJD, time, decimal year, and NEU must have the same row count.');
end
if ~isdatetime(ts.time)
    error('offset:InvalidTimeType', 'Time values must be MATLAB datetime values.');
end

mjd = ts.mjd(:);
time = ts.time(:);
decimalYear = ts.decimalYear(:);
invalidMjd = ~isfinite(mjd);
invalidTime = isnat(time);
invalidDecimalYear = ~isfinite(decimalYear);
invalidNeu = any(~isfinite(ts.neu), 2);
invalid = invalidMjd | invalidTime | invalidDecimalYear | invalidNeu;

quality = struct();
quality.originalCount = n;
quality.removedNonfiniteCount = nnz(invalid);
quality.removedOriginalIndices = find(invalid);
quality.removedOriginalLineNumbers = nan(nnz(invalid), 1);
quality.removedNonfiniteMjd = mjd(invalid);
quality.removedNonfiniteTime = time(invalid);
quality.removedNonfiniteDecimalYear = decimalYear(invalid);
quality.removedNonfiniteNeu = ts.neu(invalid, :);
quality.removalReasons = struct( ...
    'mjd', invalidMjd(invalid), ...
    'time', invalidTime(invalid), ...
    'decimalYear', invalidDecimalYear(invalid), ...
    'neu', invalidNeu(invalid));

keep = ~invalid;
ts.mjd = mjd(keep);
ts.time = time(keep);
ts.decimalYear = decimalYear(keep);
ts.neu = ts.neu(keep, :);
if isempty(ts.mjd)
    error('offset:NoFiniteObservations', ...
        'No observations remain after invalid MOM epochs are removed.');
end

[ts.mjd, order] = sort(ts.mjd);
ts.time = ts.time(order);
ts.decimalYear = ts.decimalYear(order);
ts.neu = ts.neu(order, :);

epochSpacing = diff(ts.mjd);
duplicateIndex = find(epochSpacing <= 0, 1, 'first');
if ~isempty(duplicateIndex)
    error('offset:NonIncreasingEpochs', ...
        'MOM epochs must be unique; MJD %.12g occurs more than once.', ...
        ts.mjd(duplicateIndex));
end

ts = alignMomComponents(ts);

quality.retainedCount = numel(ts.mjd);
quality.wasSorted = any(order(:) ~= (1:numel(order)).');
gapMask = epochSpacing > 1 + 1e-9;
missingDays = max(round(epochSpacing(gapMask)) - 1, 0);
gapIndices = find(gapMask);
quality.gapCount = nnz(gapMask);
quality.gapDays = sum(missingDays);
quality.gapStartMjd = ts.mjd(gapIndices);
quality.gapEndMjd = ts.mjd(gapIndices + 1);
quality.gapDurationsDays = epochSpacing(gapMask);
quality.missingDaysByGap = missingDays;

end

function ts = alignMomComponents(ts)
components = ts.metadata.mom.components;
for k = 1:numel(components)
    component = components(k);
    componentMjd = component.mjd(:);
    finiteMjd = componentMjd(isfinite(componentMjd));
    if any(diff(sort(finiteMjd)) <= 0)
        error('offset:NonIncreasingEpochs', ...
            'MOM component %d contains duplicate epochs.', k - 1);
    end

    [present, location] = ismember(ts.mjd, componentMjd);
    if any(~present)
        error('offset:MOMAlignmentFailed', ...
            'Cleaned epochs cannot be aligned with MOM component %d.', k - 1);
    end
    component.mjd = component.mjd(location);
    component.values = component.values(location);
    component.columns = component.columns(location);
    components(k) = component;
end
ts.metadata.mom.components = components;
end
