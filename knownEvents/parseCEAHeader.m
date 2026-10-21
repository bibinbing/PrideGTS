function T = parseCEAHeader(ceaFile)
%PARSECEAHEADER Extract known offset epochs from a CEA header block.
%
%   T = PARSECEAHEADER(ceaFile)
%
%   CEA (First Crust Monitoring and Application Center) GNSS coordinate
%   files annotate step epochs and amplitudes directly in the header,
%   grouped by component, e.g.:
%
%       #                    e component
%       #    offset 0.00453  +/-  0.00092 m [2011.19041]
%       #                    u component
%       #    offset -0.01006  +/-  0.00478 m [2015.90274]
%
%   Parsing stops at the first non-comment line (the data block), so the
%   (potentially large) numeric body of the file is never read into
%   memory.
%
%   Output table columns:
%     MJD       - modified Julian date (integer, whole day)
%     Component - "N" | "E" | "U"
%     Offset_m  - step amplitude (m)
%     Sigma_m   - formal uncertainty (m)
%
%   Example:
%     T = parseCEAHeader('BJFS.txt');

fid = fopen(ceaFile, 'r');
if fid < 0
    error('parseCEAHeader:fileNotFound', 'Cannot open CEA file: %s', ceaFile);
end

mjdList  = zeros(0,1);
compList = strings(0,1);
offList  = zeros(0,1);
sigList  = zeros(0,1);
curComp  = "";

compMap = containers.Map({'n','e','u'}, {'N','E','U'});

while true
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    trimmed = strtrim(tline);
    if isempty(trimmed) || trimmed(1) ~= '#'
        break;   % reached the data block: stop scanning the header
    end

    % Component switch line, e.g. "#                    e component"
    mComp = regexp(tline, '#\s*([neu])\s+component', 'tokens', 'once', 'ignorecase');
    if ~isempty(mComp)
        curComp = compMap(lower(mComp{1}));
        continue;
    end

    % Offset line, e.g. "#    offset 0.00453  +/-  0.00092 m [2011.19041]"
    mOff = regexp(tline, ...
        'offset\s+([-\d.]+)\s+\+/-\s+([-\d.]+)\s*m\s*\[\s*([\d.]+)\s*\]', ...
        'tokens', 'once');
    if ~isempty(mOff)
        offVal  = str2double(mOff{1});
        sigVal  = str2double(mOff{2});
        decYear = str2double(mOff{3});

        mjdList(end+1,1)  = decyear2mjd(decYear); %#ok<AGROW>
        compList(end+1,1) = curComp;
        offList(end+1,1)  = offVal;
        sigList(end+1,1)  = sigVal;
    end
end
fclose(fid);

MJD = mjdList;
Component = compList;
Offset_m = offList;
Sigma_m = sigList;
T = table(MJD, Component, Offset_m, Sigma_m);
if height(T) > 0
    T = sortrows(T, 'MJD');
end

end
