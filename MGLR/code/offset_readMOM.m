function ts = offset_readMOM(fileName)
% Reads a three-component MOM station group.

[folder, stem, ext] = fileparts(fileName);
if ~strcmpi(ext, '.mom')
    error('offset:MOMExtension', 'Expected a .mom file: %s.', fileName);
end

token = regexp(stem, '^(.*)_([012])$', 'tokens', 'once');
if isempty(token)
    error('offset:MOMTripletName', 'MOM input must use station_0/1/2.mom naming: %s.', fileName);
end

station = token{1};
componentFiles = cell(1, 3);
components = repmat(struct('file', '', 'headerLines', {{}}, ...
    'offsetMjd', [], 'mjd', [], 'values', [], 'columns', []), 1, 3);
for k = 1:3
    componentFiles{k} = fullfile(folder, sprintf('%s_%d%s', station, k - 1, ext));
    if ~isfile(componentFiles{k})
        error('offset:MOMMissingComponent', 'Missing MOM component file: %s.', componentFiles{k});
    end
    components(k) = readMomComponent(componentFiles{k});
end

commonMjd = components(1).mjd(:);
for k = 2:3
    commonMjd = intersect(commonMjd, components(k).mjd(:), 'stable');
end
if isempty(commonMjd)
    error('offset:MOMNoCommonEpochs', 'MOM component files for %s have no common MJD epochs.', station);
end

Y = nan(numel(commonMjd), 3);
for k = 1:3
    [isMember, loc] = ismember(commonMjd, components(k).mjd);
    if any(~isMember)
        error('offset:MOMAlignmentFailed', 'Failed to align MOM component %d for %s.', k - 1, station);
    end
    Y(:, k) = components(k).values(loc);
end

valueScaleToMeters = 1e-3;
time = offset_date(commonMjd);

ts = struct();
ts.station = upper(char(station));
ts.decimalYear = offset_decyr(time);
ts.time = time(:);
ts.mjd = commonMjd(:);
ts.neu = Y .* valueScaleToMeters;
ts.metadata = struct();
ts.metadata.mom = struct();
ts.metadata.mom.station = station;
ts.metadata.mom.componentFiles = componentFiles;
ts.metadata.mom.components = components;
ts.metadata.mom.valueScaleToMeters = valueScaleToMeters;
ts.metadata.mom.componentOrder = {'North', 'East', 'Up'};
ts.metadata.mom.knownOffsetMjd = unique(vertcat(components.offsetMjd), 'sorted');
ts.knownOffsetMjd = ts.metadata.mom.knownOffsetMjd;
ts.sourceFile = componentFiles{1};
ts.sourceFormat = 'MOM';
ts.count = numel(commonMjd);

end

function component = readMomComponent(fileName)
lines = offset_txt(fileName);
headerLines = {};
offsetMjd = [];
mjd = [];
values = [];
columns = {};

for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue;
    end
    if startsWith(line, '#')
        headerLines{end + 1, 1} = lines{i}; %#ok<AGROW>
        if ~isempty(regexpi(line, '^\s*#\s*offset(?:\s|$)', 'once'))
            token = regexpi(line, ...
                ['^\s*#\s*offset\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)', ...
                 '(?:[eEdD][+-]?\d+)?)'], 'tokens', 'once');
            if isempty(token)
                error('offset:MOMInvalidOffsetHeader', ...
                    'Invalid MOM offset header in %s at line %d: %s', ...
                    fileName, i, lines{i});
            end
            value = str2double(regexprep(token{1}, '[dD]', 'e'));
            if ~isfinite(value)
                error('offset:MOMInvalidOffsetHeader', ...
                    'MOM offset MJD must be finite in %s at line %d.', ...
                    fileName, i);
            end
            offsetMjd(end + 1, 1) = value; %#ok<AGROW>
        end
        continue;
    end

    nums = sscanf(line, '%f').';
    if numel(nums) < 2
        continue;
    end
    mjd(end + 1, 1) = nums(1); %#ok<AGROW>
    values(end + 1, 1) = nums(2); %#ok<AGROW>
    columns{end + 1, 1} = nums; %#ok<AGROW>
end

if isempty(mjd)
    error('offset:MOMEmpty', 'No MOM observations found in %s.', fileName);
end

[mjd, order] = sort(mjd);
values = values(order);
columns = columns(order);

component = struct();
component.file = fileName;
component.headerLines = headerLines;
component.offsetMjd = unique(offsetMjd, 'sorted');
component.mjd = mjd;
component.values = values;
component.columns = columns;
end
