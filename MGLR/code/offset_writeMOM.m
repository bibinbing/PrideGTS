function outFile = offset_writeMOM(ts, ~, detection, opts)
% Writes MOM components and updates their offset headers.

mom = ts.metadata.mom;
station = mom.station;

outFile = fullfile(opts.output.seriesPath, sprintf('%s_0.mom', station));
offsetMjd = [];
if ~isempty(detection) && isfield(detection, 'appliedIndices')
    offsetMjd = outputOffsetMjd(ts, detection, detection.appliedIndices);
elseif ~isempty(detection) && isfield(detection, 'acceptedIndices')
    offsetMjd = outputOffsetMjd(ts, detection, detection.acceptedIndices);
end

for k = 1:3
    component = mom.components(k);
    targetFile = fullfile(opts.output.seriesPath, sprintf('%s_%d.mom', station, k - 1));
    writeMomComponent(component, targetFile, offsetMjd);
end

end

function offsetMjd = outputOffsetMjd(ts, detection, indices)
indices = indices(:);
offsetMjd = ts.mjd(indices);
if isempty(indices) || ~isfield(detection, 'constraintMjd') || ...
        ~isfield(detection, 'acceptedIndices')
    return;
end
[present, location] = ismember(indices, detection.acceptedIndices(:));
constraintMjd = detection.constraintMjd(:);
valid = present & location >= 1 & location <= numel(constraintMjd);
replacement = nan(size(offsetMjd));
replacement(valid) = constraintMjd(location(valid));
useConstraint = valid & isfinite(replacement);
offsetMjd(useConstraint) = replacement(useConstraint);
end

function writeMomComponent(component, targetFile, offsetMjd)
fid = fopen(targetFile, 'w');
if fid < 0
    error('offset:WriteFailed', 'Cannot open %s for writing.', targetFile);
end
cleaner = onCleanup(@() fclose(fid));

headerLines = removeOffsetLines(component.headerLines);
for i = 1:numel(headerLines)
    fprintf(fid, '%s\n', headerLines{i});
end
for i = 1:numel(offsetMjd)
    fprintf(fid, '# offset %.4f\n', offsetMjd(i));
end

for i = 1:numel(component.mjd)
    values = component.columns{i};
    fprintf(fid, '%12.6f', values(1));
    for j = 2:numel(values)
        fprintf(fid, ' %.6f', values(j));
    end
    fprintf(fid, '\n');
end
end

function headerLines = removeOffsetLines(headerLines)
keep = true(numel(headerLines), 1);
for i = 1:numel(headerLines)
    keep(i) = isempty(regexpi(headerLines{i}, ...
        '^\s*#\s*offset(?:\s|$)', 'once'));
end
headerLines = headerLines(keep);
end
