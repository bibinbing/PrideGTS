function known = offset_knownEvents(ts)
% Maps MOM offset headers to retained observation indices.

known = struct('requestedMjd', zeros(0, 1), ...
    'allRequestedMjd', zeros(0, 1), ...
    'requestedMjdGroups', {cell(0, 1)}, ...
    'effectiveMjd', zeros(0, 1), 'indices', zeros(0, 1), ...
    'collapsedCount', 0);
if ~isfield(ts, 'knownOffsetMjd') || isempty(ts.knownOffsetMjd)
    return;
end
if ~isfield(ts, 'mjd') || isempty(ts.mjd)
    error('offset:MissingTimeSeriesField', ...
        'Known offset constraints require a nonempty MJD series.');
end

requestedMjd = ts.knownOffsetMjd;
if ~isnumeric(requestedMjd) || ~isreal(requestedMjd) || ...
        any(~isfinite(requestedMjd(:)))
    error('offset:InvalidKnownOffsetMjd', ...
        'Known offset MJDs must be finite real numeric values.');
end
requestedMjd = unique(requestedMjd(:), 'sorted');
mjd = ts.mjd(:);

outside = requestedMjd <= mjd(1) | requestedMjd > mjd(end);
if any(outside)
    error('offset:KnownOffsetOutsideSeries', ...
        ['Known offset MJD(s) %s do not define an identifiable step ', ...
         'strictly after the first and no later than the last observation.'], ...
        strjoin(compose('%.6f', requestedMjd(outside)), ', '));
end

indices = zeros(numel(requestedMjd), 1);
for i = 1:numel(requestedMjd)
    indices(i) = find(mjd >= requestedMjd(i), 1, 'first');
end
uniqueIndices = unique(indices, 'stable');
representativeMjd = zeros(numel(uniqueIndices), 1);
requestedMjdGroups = cell(numel(uniqueIndices), 1);
for i = 1:numel(uniqueIndices)
    groupMjd = requestedMjd(indices == uniqueIndices(i));
    requestedMjdGroups{i} = groupMjd;
    representativeMjd(i) = groupMjd(1);
end

known.requestedMjd = representativeMjd;
known.allRequestedMjd = requestedMjd;
known.requestedMjdGroups = requestedMjdGroups;
known.effectiveMjd = mjd(uniqueIndices);
known.indices = uniqueIndices;
known.collapsedCount = numel(requestedMjd) - numel(uniqueIndices);

end
