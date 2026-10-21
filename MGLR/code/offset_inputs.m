function files = offset_inputs(inputPaths)
% Collects MOM input files from files or folders.

files = {};
if isempty(inputPaths)
    return;
end
if ischar(inputPaths) || isstring(inputPaths)
    inputPaths = cellstr(inputPaths);
end

for i = 1:numel(inputPaths)
    p = char(inputPaths{i});
    if isfolder(p)
        patterns = {'*.mom', '*.MOM'};
        for j = 1:numel(patterns)
            listing = dir(fullfile(p, patterns{j}));
            for k = 1:numel(listing)
                fileName = canonicalMomPath(fullfile( ...
                    listing(k).folder, listing(k).name));
                if ~isempty(fileName)
                    files{end + 1} = fileName; %#ok<AGROW>
                end
            end
        end
    elseif isfile(p)
        fileName = canonicalMomPath(p);
        if isempty(fileName)
            warning('offset:UnsupportedInput', ...
                'Only MOM input files are supported: %s', p);
        else
            files{end + 1} = fileName; %#ok<AGROW>
        end
    else
        warning('offset:MissingInput', 'Input path does not exist: %s', p);
    end
end

files = unique(files, 'stable');

end

function fileName = canonicalMomPath(fileName)
[folder, stem, ext] = fileparts(fileName);
if ~strcmpi(ext, '.mom')
    fileName = '';
    return;
end

token = regexp(stem, '^(.*)_([012])$', 'tokens', 'once');
if isempty(token)
    return;
end

base = token{1};
candidate = fullfile(folder, [base '_0' ext]);
if isfile(candidate)
    fileName = candidate;
end
end
