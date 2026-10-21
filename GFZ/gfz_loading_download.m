function report = gfz_loading_download(cfg, parent)
%GFZ_LOADING_DOWNLOAD GFZ CF point series for App Designer (R2020b+ syntax).
% Required cfg: StationFile, SaveRoot, StartDate, EndDate, Types, Directions.
% EndDate is EXCLUSIVE. Missing station dates use the corresponding default.
% Optional: UseProcessPool=true, ChunkDays=90, Timeout=60, MaxAttempts=2.
% UI objects remain on the client. Only HTTP is dispatched to a process pool.
% Output: NetCDF/merged/STATION_TYPE.nc; incomplete runs get separate files.
if nargin < 2, parent = []; end
cfg = defaults_(cfg);
report = struct('Canceled', false, 'MergedFiles', {{}}, ...
    'PartialFiles', {{}}, 'Failures', {{}}, 'Requests', 0, ...
    'CacheHits', 0, 'LogFile', '');
dialog = []; pool = []; activeFuture = []; logId = -1;
networkFailures = 0; circuitOpen = false; done = 0; total = 1; budget = 0;
resources = containers.Map('KeyType','char','ValueType','any');
resources('dialog') = []; resources('future') = []; resources('log') = -1;
cleaner = onCleanup(@() cleanup_resources_(resources));
stations = stations_(cfg);
root = fullfile(char(cfg.SaveRoot), 'NetCDF');
mkdir_checked_(root);
mkdir_checked_(fullfile(root, 'merged'));
mkdir_checked_(fullfile(root, 'logs'));
runTag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.LogFile = fullfile(root, 'logs', ['GFZ_' runTag '.log']);
logId = fopen(report.LogFile, 'w', 'n', 'UTF-8');
if logId < 0, error('GFZ:IO', 'Cannot create log: %s', report.LogFile); end
resources('log') = logId;
if ~isempty(parent)
    dialog = uiprogressdlg(parent, 'Title', 'GFZ loading download', ...
        'Message', 'Reading GFZ annual catalogs...', 'Cancelable', 'on');
    resources('dialog') = dialog;
end
log_('Start. Dates are UTC calendar dates; end is exclusive.');
% Create at most one new worker. Reuse an existing PROCESS pool, never close it.
if cfg.UseProcessPool && license('test', 'Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isa(pool, 'parallel.ThreadPool')
            pool = [];
            log_('Existing thread pool: using synchronous HTTP.');
        elseif isempty(pool)
            ui_('Starting one process worker for responsive downloads...');
            pool = parpool('local', 1);
        end
    catch ME
        pool = [];
        log_(['Process pool unavailable; synchronous HTTP: ' ME.message]);
    end
else
    log_('Synchronous HTTP; cancellation is checked between requests.');
end
catalogs = containers.Map('KeyType', 'char', 'ValueType', 'any');
total = height(stations) * numel(cfg.Types);
try
    for s = 1:height(stations)
        for t = 1:numel(cfg.Types)
            check_cancel_();
            kind = cfg.Types{t}; name = char(stations.Name(s));
            label = [name ' ' kind];
            files = {}; issuesBefore = numel(report.Failures);
            firstYear = year(stations.Start(s));
            lastYear = year(stations.End(s) - seconds(1));
            try
                for yr = firstYear:lastYear
                    check_cancel_();
                    if yr < 1976
                        fail_([label ': requested year before configured 1976 coverage.']);
                        continue
                    end
                    period = 1 + (yr >= 2000) + (yr >= 2010);
                    mapKey = sprintf('%s_%d', kind, period);
                    if ~isKey(catalogs, mapKey)
                        id = entry_id_(kind, period);
                        url = ['https://esmdata.gfz.de/repository/entry/show?entryid=' id];
                        ui_(['Reading GFZ catalog: ' kind]);
                        try
                            response = get_(url);
                            text = native2unicode(response.Bytes(:)', 'UTF-8');
                            entries = catalog_(text, kind);
                            catalogs(mapKey) = entries;
                        catch ME
                            if strcmp(ME.identifier, 'GFZ:Canceled'), rethrow(ME); end
                            catalogs(mapKey) = struct([]);
                            fail_([label ': catalog unavailable: ' ME.message]);
                        end
                    end
                    entries = catalogs(mapKey);
                    if isempty(entries)
                        fail_(sprintf('%s: annual catalog unavailable for %d.', label, yr));
                        continue
                    end
                    idx = find([entries.Year] == yr);
                    if numel(idx) ~= 1
                        fail_(sprintf('%s: no unique annual CF entry for %d.', label, yr));
                        continue
                    end
                    entry = entries(idx);
                    a = max(posixtime(stations.Start(s)), utc_seconds_(sprintf('%04d-01-01', yr)));
                    b = min(posixtime(stations.End(s)), utc_seconds_(sprintf('%04d-01-01', yr + 1)));
                    % Use the source time grid, including HYDL's 12:00 phase.
                    wantedFirst = entry.First + ceil((a-entry.First)/entry.Step - 1e-9)*entry.Step;
                    wantedLast = entry.First + (ceil((b-entry.First)/entry.Step - 1e-9)-1)*entry.Step;
                    if wantedFirst > wantedLast, continue; end
                    if wantedFirst < entry.First || wantedLast > entry.Last
                        fail_(sprintf('%s %d: requested coverage exceeds published %s to %s.', ...
                            label, yr, iso_(entry.First), iso_(entry.Last)));
                    end
                    a = max(wantedFirst, entry.First); b = min(wantedLast, entry.Last);
                    nPerChunk = max(1, floor(cfg.ChunkDays*86400/entry.Step));
                    while a <= b
                        z = min(b, a + (nPerChunk-1)*entry.Step);
                        task = struct('Name', name, 'Kind', kind, 'Entry', entry, ...
                            'Lat', stations.Lat(s), 'Lon', stations.Lon(s), ...
                            'First', a, 'Last', z, 'Directions', {cfg.Directions});
                        budget = 7; % Shared across one top-level chunk and all children.
                        got = chunk_(task, 0);
                        files = [files, got]; %#ok<AGROW>
                        a = z + entry.Step;
                    end
                end
            catch ME
                if strcmp(ME.identifier, 'GFZ:Canceled')
                    report.Canceled = true;
                else
                    fail_([label ': ' ME.message]);
                end
            end
            % Cancellation keeps completed chunks and merges them as PARTIAL.
            isPartial = report.Canceled || numel(report.Failures) > issuesBefore;
            if ~isempty(files)
                ui_(['Merging ' label '...']);
                if isPartial
                    out = fullfile(root, 'merged', [name '_' kind '_partial_' runTag '.nc']);
                else
                    out = fullfile(root, 'merged', [name '_' kind '.nc']);
                end
                try
                    merge_(files, out, cfg.Directions, kind, stations.Start(s), ...
                        stations.End(s), isPartial, name);
                    if isPartial
                        report.PartialFiles{end+1} = out;
                    else
                        report.MergedFiles{end+1} = out;
                    end
                    log_(['Saved: ' out]);
                catch ME
                    fail_([label ': merge failed: ' ME.message]);
                end
            elseif ~report.Canceled
                fail_([label ': no validated data available to merge.']);
            end
            done = done + 1;
            if report.Canceled, break; end
        end
        if report.Canceled, break; end
    end
catch ME
    if strcmp(ME.identifier, 'GFZ:Canceled')
        report.Canceled = true;
    else
        fail_(ME.message);
    end
end
log_(sprintf('Finished: %d complete, %d partial, %d issues, %d requests, %d cache hits.', ...
    numel(report.MergedFiles), numel(report.PartialFiles), ...
    numel(report.Failures), report.Requests, report.CacheHits));

    function files = chunk_(task, depth)
        check_cancel_();
        files = {};
        identity = sprintf('%s|%s|%.17g|%.17g|%s', task.Kind, ...
            task.Entry.Id, task.Lat, task.Lon, strjoin(task.Directions, ','));
        cacheDir = fullfile(root, 'parts_v2', digest_(identity));
        mkdir_checked_(cacheDir);
        file = fullfile(cacheDir, sprintf('%.0f_%.0f.nc', task.First, task.Last));
        if isfile(file)
            try
                validate_(file, task, true);
                report.CacheHits = report.CacheHits + 1;
                files = {file}; return
            catch ME
                log_(['Cache rejected: ' file ': ' ME.message]);
            end
        end
        % Parents that previously split should reuse children immediately.
        splitMarker = [file '.split'];
        if isfile(splitMarker) && can_split_(task, depth)
            files = children_(task, depth); return
        end
        ui_(sprintf('%s %s: %s to %s', task.Name, task.Kind, ...
            iso_(task.First), iso_(task.Last)));
        url = point_url_(task);
        tmp = [tempname(cacheDir) '.nc'];
        tmpCleaner = onCleanup(@() delete_if_(tmp));
        try
            response = get_(url);
            write_bytes_(tmp, response.Bytes);
            validate_(tmp, task, false);
            ncwriteatt(tmp, '/', 'gfz_entry_id', task.Entry.Id);
            ncwriteatt(tmp, '/', 'gfz_requested_latitude', task.Lat);
            ncwriteatt(tmp, '/', 'gfz_requested_longitude', task.Lon);
            ncwriteatt(tmp, '/', 'gfz_request_url', url);
            ncwriteatt(tmp, '/', 'gfz_source_name', task.Entry.Name);
            ncwriteatt(tmp, '/', 'gfz_unit_policy', 'Use explicit metres or GFZ displacement product metre convention when omitted');
            move_checked_(tmp, file);
            files = {file};
        catch ME
            if strcmp(ME.identifier, 'GFZ:Canceled'), rethrow(ME); end
            % Only request size/timeouts justify splitting. Never split HTML,
            % time-axis errors, bad credentials, 404, rate limiting, or 500.
            if any(strcmp(ME.identifier, {'GFZ:Timeout','GFZ:HTTP413','GFZ:HTTP504'})) ...
                    && can_split_(task, depth) && budget > 0 && ~circuitOpen
                log_(['Splitting once smaller: ' ME.message]);
                write_bytes_(splitMarker, uint8('split'));
                files = children_(task, depth);
            else
                fail_(sprintf('%s %s %s to %s: %s', task.Name, task.Kind, ...
                    iso_(task.First), iso_(task.Last), ME.message));
            end
        end
    end

    function files = children_(task, depth)
        n = round((task.Last-task.First)/task.Entry.Step) + 1;
        left = task; right = task;
        left.Last = task.First + (floor(n/2)-1)*task.Entry.Step;
        right.First = left.Last + task.Entry.Step;
        files = chunk_(left, depth+1);
        % Preserve left results if cancellation occurs during the right half.
        try
            files = [files, chunk_(right, depth+1)];
        catch ME
            if strcmp(ME.identifier, 'GFZ:Canceled')
                report.Canceled = true;
            else
                rethrow(ME)
            end
        end
    end

    function tf = can_split_(task, depth)
        n = round((task.Last-task.First)/task.Entry.Step) + 1;
        tf = depth < 2 && floor(n/2)*task.Entry.Step >= 31*86400;
    end

    function response = get_(url)
        if circuitOpen, error('GFZ:Circuit', 'Network circuit open; see earlier log entry.'); end
        for attempt = 1:cfg.MaxAttempts
            check_cancel_();
            % budget exists only while downloading a station chunk.
            isPoint = contains(url, 'output=data.gridaspoint');
            if isPoint
                if budget <= 0, error('GFZ:Budget', 'Chunk request budget exhausted.'); end
                budget = budget - 1;
            end
            report.Requests = report.Requests + 1;
            log_(sprintf('Request %d attempt %d: %s', report.Requests, attempt, url));
            delay = 2^(attempt-1) + rand;
            try
                if isempty(pool)
                    response = cfg.HTTPFunction(url, cfg.Timeout);
                else
                    activeFuture = parfeval(pool, cfg.HTTPFunction, 1, url, cfg.Timeout);
                    resources('future') = activeFuture;
                    while ~strcmp(activeFuture.State, 'finished')
                        check_cancel_(); pause(0.1);
                    end
                    response = fetchOutputs(activeFuture);
                    delete(activeFuture); activeFuture = []; resources('future') = [];
                end
            catch ME
                if ~isempty(activeFuture)
                    cancel(activeFuture); delete(activeFuture); activeFuture = []; resources('future') = [];
                end
                if strcmp(ME.identifier, 'GFZ:Canceled'), rethrow(ME); end
                message = lower([ME.identifier ' ' ME.message]);
                isTimeout = contains(message, 'timeout') || contains(message, 'timed out');
                retryable = isTimeout || contains(message, 'connection') || contains(message, 'connect');
                networkFailures = networkFailures + 1;
                if networkFailures >= 3, circuitOpen = true; end
                if ~retryable || attempt == cfg.MaxAttempts || circuitOpen
                    if isTimeout, error('GFZ:Timeout', '%s', ME.message); end
                    rethrow(ME)
                end
                log_(['Transport retry: ' ME.message]);
                wait_(delay); continue
            end
            status = response.Status;
            if status == 200
                networkFailures = 0; return
            end
            log_(sprintf('HTTP %d', status));
            if any(status == [401 403]), circuitOpen = true; end
            transient = any(status == [408 429 500 502 503 504]);
            if transient
                networkFailures = networkFailures + 1;
                if networkFailures >= 3, circuitOpen = true; end
            end
            delay = max(delay, response.RetryAfter);
            if ~transient || attempt == cfg.MaxAttempts || circuitOpen || delay > 60
                % Stop session-wide for exhausted rate limiting/service outages.
                if any(status == [429 503]) || delay > 60, circuitOpen = true; end
                error(sprintf('GFZ:HTTP%d', status), ...
                    'HTTP %d; retry-after %.1f s. URL: %s', status, response.RetryAfter, url);
            end
            wait_(delay);
        end
    end

    function wait_(duration)
        start = tic;
        while toc(start) < duration
            check_cancel_(); pause(min(0.1, max(0, duration-toc(start))));
        end
    end
    function check_cancel_()
        drawnow limitrate;
        if report.Canceled || (~isempty(parent) && ~isvalid(parent)) || ...
                (~isempty(dialog) && (~isvalid(dialog) || dialog.CancelRequested))
            error('GFZ:Canceled', 'Canceled by user.');
        end
    end
    function ui_(message)
        if ~isempty(dialog) && isvalid(dialog)
            dialog.Message = message; dialog.Value = min(1, done/total);
        end
        drawnow limitrate;
    end
    function log_(message)
        if logId >= 0
            fprintf(logId, '%s %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), message);
        end
    end
    function fail_(message)
        report.Failures{end+1} = message; log_(['ISSUE: ' message]);
    end
end

function key = digest_(text)
% Compact identity avoids Windows/NetCDF failures with deeply nested paths.
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(unicode2native(text,'UTF-8'));
bytes = typecast(md.digest(),'uint8');
key = lower(reshape(dec2hex(bytes,2).',1,[])); key = key(1:32);
end

function cleanup_resources_(resources)
f = resources('future');
if ~isempty(f) && isvalid(f), cancel(f); delete(f); end
d = resources('dialog');
if ~isempty(d) && isvalid(d), close(d); end
fid = resources('log');
if fid >= 0, fclose(fid); end
end

function cfg = defaults_(cfg)
cfg.SaveRoot = text_path_(cfg.SaveRoot, 'Please select a download save directory first.');
d = struct('UseProcessPool', true, 'ChunkDays', 90, 'Timeout', 60, ...
    'MaxAttempts', 2, 'HTTPFunction', @gfz_http_get);
names = fieldnames(d);
for k = 1:numel(names)
    if ~isfield(cfg, names{k}), cfg.(names{k}) = d.(names{k}); end
end
cfg.Types = unique(cellstr(string(cfg.Types)), 'stable');
cfg.Directions = sort(unique(cellstr(string(cfg.Directions))));
assert(~isempty(cfg.Types) && all(ismember(cfg.Types, {'HYDL','NTAL','NTOL'})), ...
    'GFZ:Input', 'Select HYDL, NTAL and/or NTOL.');
assert(~isempty(cfg.Directions) && all(ismember(cfg.Directions, {'duEW','duNS','duV'})), ...
    'GFZ:Input', 'Select at least one displacement direction.');
validateattributes(cfg.ChunkDays, {'numeric'}, {'scalar','finite','positive','<=',366});
validateattributes(cfg.Timeout, {'numeric'}, {'scalar','finite','positive'});
validateattributes(cfg.MaxAttempts, {'numeric'}, {'scalar','integer','>=',1,'<=',3});
end

function path = text_path_(value, message)
assert(ischar(value) || isstring(value), 'GFZ:Input', '%s', message);
path = strtrim(string(value));
assert(isscalar(path) && ~ismissing(path) && strlength(path)>0, ...
    'GFZ:Input', '%s', message);
path = char(path);
end

function st = stations_(cfg)
cfg.StationFile = text_path_(cfg.StationFile, 'Please select a station information file first.');
assert(isfile(cfg.StationFile), 'GFZ:Input', 'Station file does not exist.');
opts = detectImportOptions(cfg.StationFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'MultipleDelimsAsOne', false, 'ReadVariableNames', false);
opts = setvartype(opts, opts.VariableNames, 'string');
raw = readtable(cfg.StationFile, opts);
assert(height(raw)>0 && width(raw)>=3, 'GFZ:Input', ...
    'Use a headerless tab-separated station file: name, latitude, longitude, [start], [end].');
n = height(raw); name = strtrim(string(raw{:,1}));
lat = str2double(string(raw{:,2})); lon = str2double(string(raw{:,3}));
assert(all(isfinite(lat) & abs(lat)<=90) && all(isfinite(lon) & lon>=-180 & lon<=360), ...
    'GFZ:Input', 'Latitude must be in [-90,90], longitude in [-180,360].');
lon = mod(lon+180,360)-180;
assert(~any(ismissing(name)), 'GFZ:Input', 'Station names cannot be empty.');
for k = 1:n
    assert(~isempty(regexp(char(name(k)), '^[A-Za-z0-9_-]+$', 'once')) && ...
        isempty(regexpi(char(name(k)), '^(CON|PRN|AUX|NUL|COM[0-9]|LPT[0-9])$', 'once')), ...
        'GFZ:Input', 'Use safe station names containing letters, digits, _ or -.');
end
assert(numel(unique(lower(name)))==n, 'GFZ:Input', 'Station names must be unique (case insensitive).');
a = NaT(n,1,'TimeZone','UTC'); b = a;
for k = 1:n
    x = ""; y = "";
    if width(raw)>=4, x = strtrim(string(raw{k,4})); end
    if width(raw)>=5, y = strtrim(string(raw{k,5})); end
    a(k) = date_(x, cfg.StartDate); b(k) = date_(y, cfg.EndDate);
end
assert(all(a<b), 'GFZ:Input', 'End date must be later than start date; end is exclusive.');
st = table(name,lat,lon,a,b,'VariableNames',{'Name','Lat','Lon','Start','End'});
end

function d = date_(value, fallback)
if ismissing(value) || strlength(value)==0
    assert(isdatetime(fallback) && isscalar(fallback) && ~isnat(fallback), ...
        'GFZ:Input', 'Select the missing default date.');
    value = string(fallback, 'yyyy-MM-dd');
end
assert(~isempty(regexp(char(value), '^\d{4}-\d{2}-\d{2}$', 'once')), ...
    'GFZ:Input', 'Station dates must have format yyyy-MM-dd.');
d = datetime(value, 'InputFormat','yyyy-MM-dd','TimeZone','UTC');
assert(~isnat(d) && string(d,'yyyy-MM-dd')==value, 'GFZ:Input', 'Invalid calendar date.');
end

function id = entry_id_(kind, period)
switch kind
    case 'HYDL'
        ids = {'852a85f3-01e7-40fb-9afe-b2cad5518028', ...
            'f467021f-2f69-4624-9e09-a100961aba09','d5fbcb93-a44e-4b0c-aa0d-6db056597ebe'};
    case 'NTAL'
        ids = {'17551ce9-ae62-45f9-8f95-584b99fb3ad9', ...
            'af9d84c4-d6be-4342-ad9a-12e87c909880','538d0cd7-8ef0-4670-b85f-4da6f6f2d6d1'};
    case 'NTOL'
        ids = {'be9fcbef-684f-4844-95b2-b08900d0dae1', ...
            'd044c8e5-751f-42f7-ab85-254e0f811d31','43e64376-0147-4b26-83c2-5204659a684c'};
end
id = ids{period};
end

function entries = catalog_(html, kind)
% Parse the actual RAMADDA embedded JSON table, NOT arbitrary HTML dates.
token = regexp(html, 'var\s+entries_\w+\s*=\s*(\[[\s\S]*?\])\s*;', 'tokens', 'once');
assert(~isempty(token), 'GFZ:Catalog', 'GFZ annual catalog format changed; inspect the URL in the log.');
rows = jsondecode(token{1});
entries = struct('Id',{},'Name',{},'Year',{},'First',{},'Last',{},'Step',{});
for k = 1:numel(rows)
    if iscell(rows), row = rows{k}; else, row = rows(k); end
    if ~isfield(row,'name'), continue; end
    match = regexp(row.name, ['^ESMGFZ_' kind '_cf_[\w.]+_(\d+)h_(\d{4})\.nc$'], 'tokens','once');
    if isempty(match), continue; end
    step = str2double(match{1})*3600;
    a = utc_seconds_(row.startDate); b = utc_seconds_(row.endDate);
    assert(isfinite(step) && step>0 && a<=b, 'GFZ:Catalog', 'Invalid annual time metadata.');
    entries(end+1) = struct('Id',row.id,'Name',row.name,'Year',str2double(match{2}), ...
        'First',a,'Last',b,'Step',step); %#ok<AGROW>
end
assert(~isempty(entries), 'GFZ:Catalog', 'No supported annual CF files found in catalog.');
end

function url = point_url_(task)
url = sprintf(['https://esmdata.gfz.de/repository/entry/show?submit=Get%%20Data' ...
    '&output=data.gridaspoint&entryid=%s&location.latitude=%.10f&location.longitude=%.10f' ...
    '&calendar=proleptic_gregorian&fromdate_exact=%s&todate_exact=%s&format=netcdf'], ...
    task.Entry.Id, task.Lat, task.Lon, encode_time_(task.First), encode_time_(task.Last));
for k = 1:numel(task.Directions), url = [url '&var=' task.Directions{k}]; end %#ok<AGROW>
end

function s = encode_time_(t)
s = strrep(iso_(t), 'T', '%20'); s = strrep(s, ':', '%3A');
end
function s = iso_(t)
s = char(datetime(t,'ConvertFrom','posixtime','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss'));
end

function t = utc_seconds_(text)
% CF epoch: date alone, T/space, fractional seconds, Z/UTC or numeric offset.
text = strtrim(char(text));
pat = ['^(\d{4})-(\d{1,2})-(\d{1,2})' ...
    '(?:[T ](\d{1,2}):(\d{1,2})(?::(\d+(?:\.\d+)?))?)?' ...
    '\s*(Z|UTC|[+-]\d{2}:?\d{2})?$'];
assert(~isempty(regexp(text,pat,'once')), 'GFZ:Time', 'Unsupported time origin: %s', text);
% MATLAB regexp omits nested optional tokens on some releases. Extract
% the timezone and numerical fields independently after validating syntax.
zone = regexp(text,'(Z|UTC|[+-]\d{2}:?\d{2})$','match','once');
if ~isempty(zone), text = text(1:end-numel(zone)); end
m = regexp(text,'\d+(?:\.\d+)?','match');
v = zeros(1,6);
for k = 1:min(6,numel(m))
    if ~isempty(m{k}), v(k) = str2double(m{k}); end
end
d = datetime(v(1),v(2),v(3),v(4),v(5),v(6),'TimeZone','UTC');
t = posixtime(d);
if ~isempty(zone) && any(zone(1)==['+' '-'])
    offset = strrep(zone, ':', '');
    z = str2double(offset(2:3))*3600 + str2double(offset(4:5))*60;
    if offset(1)=='-', z = -z; end
    t = t-z;
end
assert(isfinite(t), 'GFZ:Time', 'Invalid time origin.');
end

function series = read_series_(file, directions)
info = ncinfo(file); names = {info.Variables.Name};
assert(all(ismember([{'time','latitude','longitude'}, directions],names)), ...
    'GFZ:Schema', 'Required time/coordinate/displacement variable is missing.');
units = char(ncreadatt(file,'time','units'));
m = regexpi(strtrim(units),'^(seconds?|minutes?|hours?|days?)\s+since\s+(.+)$','tokens','once');
assert(~isempty(m), 'GFZ:Time', 'Unsupported CF time units: %s', units);
scale = 1;
switch lower(m{1}(1))
    case 'm', scale = 60;
    case 'h', scale = 3600;
    case 'd', scale = 86400;
end
timeInfo = info.Variables(strcmp(names,'time'));
attrs = {timeInfo.Attributes.Name};
if any(strcmp(attrs,'calendar'))
    calendar = char(ncreadatt(file,'time','calendar'));
    assert(any(strcmp(calendar,{'standard','gregorian','proleptic_gregorian'})), ...
        'GFZ:Time', 'Unsupported calendar: %s', calendar);
end
series.Time = utc_seconds_(m{2}) + double(ncread(file,'time'))*scale;
series.Time = series.Time(:);
series.Lat = double(ncread(file,'latitude')); series.Lon = double(ncread(file,'longitude'));
assert(isscalar(series.Lat) && isfinite(series.Lat) && ...
    isscalar(series.Lon) && isfinite(series.Lon), 'GFZ:Schema', 'Expected one finite station coordinate.');
series.Values = zeros(numel(series.Time),numel(directions));
series.Units = cell(1,numel(directions));
for k = 1:numel(directions)
    v = info.Variables(strcmp(names,directions{k}));
    x = double(ncread(file,directions{k}));
    assert(numel(x)==numel(series.Time) && all(isfinite(x(:))), ...
        'GFZ:Data', 'Length mismatch or missing displacement values in %s.', directions{k});
    assert(~any(abs(x(:))>1e30), 'GFZ:Data', 'Unlabelled GFZ fill value in %s.', directions{k});
    for j = 1:numel(v.Attributes)
        if any(strcmp(v.Attributes(j).Name, {'missing_value','_FillValue'}))
            assert(~any(ismember(x(:),double(v.Attributes(j).Value))), ...
                'GFZ:Data', 'Fill values found in %s.', directions{k});
        end
    end
    namesA = {v.Attributes.Name};
    if any(strcmp(namesA,'units'))
        unit = char(ncreadatt(file,directions{k},'units'));
    elseif any(strcmp(namesA,'long_name')) && contains(char(ncreadatt(file,directions{k},'long_name')),'[m]')
        unit = 'm'; % Verified GFZ point files put [m] in long_name.
    else
        % NTAL/NTOL annual point samples omit both units and descriptive
        % long_name. Retain the GFZ displacement-product metre convention.
        % This assumption is recorded in cache and merged-file metadata.
        unit = 'm';
    end
    assert(strcmp(strtrim(unit),'m'), 'GFZ:Units', 'Unexpected displacement unit: %s', unit);
    series.Values(:,k) = x(:); series.Units{k} = 'm';
end
end

function validate_(file, task, cache)
s = read_series_(file, task.Directions);
expected = (task.First:task.Entry.Step:task.Last)';
assert(numel(s.Time)==numel(expected) && all(abs(s.Time-expected)<0.01), ...
    'GFZ:Coverage', 'Time axis does not match requested samples (gap, duplicate or boundary mismatch).');
assert(abs(s.Lat-task.Lat)<=0.51 && abs(mod(s.Lon-task.Lon+180,360)-180)<=0.51, ...
    'GFZ:Coordinate', 'Returned location differs from requested station.');
if cache
    assert(strcmp(char(ncreadatt(file,'/','gfz_entry_id')),task.Entry.Id) && ...
        abs(double(ncreadatt(file,'/','gfz_requested_latitude'))-task.Lat)<1e-10 && ...
        abs(double(ncreadatt(file,'/','gfz_requested_longitude'))-task.Lon)<1e-10, ...
        'GFZ:Cache', 'Cache request identity mismatch.');
end
end

function merge_(files, output, directions, kind, startDate, endDate, partial, stationName)
files = unique(files,'stable'); times = cell(numel(files),1); values = times;
sources = cell(numel(files),1); location = [];
for k = 1:numel(files)
    s = read_series_(files{k},directions);
    if isempty(location), location = [s.Lat,s.Lon]; end
    assert(max(abs(location-[s.Lat,s.Lon]))<1e-8, 'GFZ:Coordinate', 'Coordinates differ across chunks.');
    times{k} = s.Time; values{k} = s.Values;
    sources{k} = char(ncreadatt(files{k},'/','gfz_source_name'));
end
t = vertcat(times{:}); x = vertcat(values{:});
[t, order] = sort(t); x = x(order,:);
duplicates = find(diff(t)==0);
for k = duplicates(:)'
    assert(all(abs(x(k,:)-x(k+1,:))<1e-10), 'GFZ:Conflict', 'Conflicting values at duplicate timestamp.');
end
[t,keep] = unique(t,'stable'); x = x(keep,:);
inside = t>=posixtime(startDate) & t<posixtime(endDate);
t = t(inside); x = x(inside,:);
assert(~isempty(t), 'GFZ:Data', 'No data in the requested interval.');
tmp = [tempname(fileparts(output)) '.nc'];
cleaner = onCleanup(@() delete_if_(tmp));
n = numel(t);
nccreate(tmp,'time','Dimensions',{'obs',n},'Datatype','double','Format','netcdf4_classic');
% Keep the original application's hours/days time-unit convention.
if strcmp(kind,'HYDL')
    timeUnits = 'days since 1976-1-1 12:00:00';
    numericTime = (t-utc_seconds_('1976-01-01T12:00:00Z'))/86400;
else
    timeUnits = 'hours since 1976-1-1 00:00:00';
    numericTime = (t-utc_seconds_('1976-01-01T00:00:00Z'))/3600;
end
ncwrite(tmp,'time',numericTime);
ncwriteatt(tmp,'time','units',timeUnits);
ncwriteatt(tmp,'time','calendar','proleptic_gregorian');
ncwriteatt(tmp,'time','standard_name','time');
for k = 1:2
    names = {'latitude','longitude'}; units = {'degrees_north','degrees_east'};
    nccreate(tmp,names{k},'Dimensions',{'station',1},'Datatype','double');
    ncwrite(tmp,names{k},location(k)); ncwriteatt(tmp,names{k},'units',units{k});
end
nccreate(tmp,'station_id','Dimensions',{'station_name_strlen',numel(stationName),'station',1},'Datatype','char');
ncwrite(tmp,'station_id',stationName(:));
ncwriteatt(tmp,'station_id','cf_role','timeseries_id');
nccreate(tmp,'stationIndex','Dimensions',{'obs',n},'Datatype','int32');
ncwrite(tmp,'stationIndex',zeros(n,1,'int32'));
ncwriteatt(tmp,'stationIndex','instance_dimension','station');
for k = 1:numel(directions)
    nccreate(tmp,directions{k},'Dimensions',{'obs',n},'Datatype','double','DeflateLevel',1);
    ncwrite(tmp,directions{k},x(:,k)); ncwriteatt(tmp,directions{k},'units','m');
    ncwriteatt(tmp,directions{k},'coordinates','latitude longitude time');
end
ncwriteatt(tmp,'/','Conventions','CF-1.6');
ncwriteatt(tmp,'/','featureType','timeSeries');
ncwriteatt(tmp,'/','loading_type',kind);
ncwriteatt(tmp,'/','reference_frame','CF');
ncwriteatt(tmp,'/','source_files',strjoin(unique(sources,'stable'),'; '));
ncwriteatt(tmp,'/','spatial_method','GFZ point service: nearest neighbour');
ncwriteatt(tmp,'/','unit_policy','Use explicit metres or GFZ displacement product metre convention when omitted');
ncwriteatt(tmp,'/','requested_start',char(string(startDate,'yyyy-MM-dd')));
ncwriteatt(tmp,'/','requested_end_exclusive',char(string(endDate,'yyyy-MM-dd')));
ncwriteatt(tmp,'/','partial',int32(partial));
ncwriteatt(tmp,'/','time_coverage_start',[iso_(t(1)) 'Z']);
ncwriteatt(tmp,'/','time_coverage_end',[iso_(t(end)) 'Z']);
% Verify the temporary result before replacing any previous successful result.
check = read_series_(tmp,directions);
assert(isequal(check.Time,t) && isequal(check.Values,x), 'GFZ:IO', 'Merged read-back mismatch.');
move_checked_(tmp,output);
end

function write_bytes_(file, bytes)
fid = fopen(file,'wb');
assert(fid>=0, 'GFZ:IO', 'Cannot write %s.',file);
guard = onCleanup(@() fclose(fid));
n = fwrite(fid,bytes,'uint8');
assert(n==numel(bytes), 'GFZ:IO', 'Incomplete disk write: %s.', file);
end
function mkdir_checked_(path)
if ~isfolder(path)
    [ok,message] = mkdir(path); assert(ok,'GFZ:IO','%s',message);
end
end
function move_checked_(source,target)
[ok,message] = movefile(source,target,'f'); assert(ok,'GFZ:IO','%s',message);
end
function delete_if_(file)
if isfile(file), delete(file); end
end
