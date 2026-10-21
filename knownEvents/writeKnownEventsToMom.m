function writeKnownEventsToMom(momFile, eventsTable, ~, outFile)
%WRITEKNOWNEVENTSTOMOM Write standard offset comments to a .mom file.
%
%   writeKnownEventsToMom(momFile, eventsTable, compLabel, outFile)
%
%   Prepends one standard MOM comment per known event:
%
%       # offset <MJD>
%
%   Only events inside the MOM file's own MJD span are written. Duplicate
%   epochs are collapsed. Re-running the function replaces existing
%   "# offset" lines and removes legacy PrideGTS KNOWN_OFFSET blocks.
%
%   eventsTable only needs a numeric MJD column. compLabel is retained in
%   the function signature for compatibility but is not written.

BEGIN_TAG = '# ===== PrideGTS KNOWN_OFFSET BEGIN =====';
END_TAG   = '# ===== PrideGTS KNOWN_OFFSET END =====';

raw = readmatrix(momFile, 'FileType', 'text', 'CommentStyle', '#');
if isempty(raw)
    error('writeKnownEventsToMom:emptyFile', 'No numeric data found in: %s', momFile);
end
mjdValues = raw(isfinite(raw(:, 1)), 1);
if isempty(mjdValues)
    error('writeKnownEventsToMom:emptyFile', ...
        'No finite MJD values found in: %s', momFile);
end
mjdMin = min(mjdValues);
mjdMax = max(mjdValues);

if ~isempty(eventsTable) && height(eventsTable) > 0
    keep = eventsTable.MJD >= mjdMin & eventsTable.MJD <= mjdMax;
    eventsTable = eventsTable(keep, :);
end
offsets = unique(eventsTable.MJD, 'sorted');

% Read the source text before opening the output, which also makes an
% in-place call (momFile == outFile) safe from immediate truncation.
fidIn = fopen(momFile, 'r');
if fidIn < 0
    error('writeKnownEventsToMom:cannotRead', ...
        'Cannot open input file: %s', momFile);
end
rawLines = textscan(fidIn, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fidIn);
lines = rawLines{1};

fidOut = fopen(outFile, 'w');
if fidOut < 0
    error('writeKnownEventsToMom:cannotWrite', 'Cannot open output file: %s', outFile);
end

for i = 1:numel(offsets)
    fprintf(fidOut, '# offset %.6f\n', offsets(i));
end

% Copy the original content while removing both the previous PrideGTS block
% and any standalone offset comments, so repeated writes stay idempotent.
inOldBlock = false;
for k = 1:numel(lines)
    trimmed = strtrim(lines{k});
    if strcmp(trimmed, BEGIN_TAG)
        inOldBlock = true;
    elseif strcmp(trimmed, END_TAG)
        inOldBlock = false;
    elseif ~inOldBlock && isempty(regexp(trimmed, ...
            '^#\s*offset\s+[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?:\s.*)?$', ...
            'once', 'ignorecase'))
        fprintf(fidOut, '%s\n', lines{k});
    end
end
fclose(fidOut);

end
