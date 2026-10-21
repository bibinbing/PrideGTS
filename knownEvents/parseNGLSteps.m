function T = parseNGLSteps(stepsFile, station)
%PARSENGLSTEPS Extract known offset epochs for one station from steps.txt.
%
%   T = PARSENGLSTEPS(stepsFile, station)
%
%   Reads an NGL-format steps catalogue
%   (https://geodesy.unr.edu/NGLStationPages/steps.txt) and returns every
%   entry for the given 4-character station code, converted to MJD.
%
%   Each line of steps.txt has one of two shapes:
%     <site> <yymmmdd> 1 <Equipment_Change_Description>
%     <site> <yymmmdd> 2 <thresholdKm> <distanceKm> <magnitude> <eventID>
%
%   Code 1 = antenna/receiver/radome/elevation-mask metadata change.
%   Code 2 = candidate co-seismic offset (only listed when the station
%            lies within NGL's magnitude-dependent distance threshold).
%
%   Output table columns:
%     MJD         - modified Julian date (integer, whole day)
%     Type        - "EQUIP" | "EQ"
%     Description - free-text description of the event
%
%   Example:
%     T = parseNGLSteps('steps.txt', 'HKSC');

station = upper(strtrim(station));

fid = fopen(stepsFile, 'r');
if fid < 0
    error('parseNGLSteps:fileNotFound', 'Cannot open steps file: %s', stepsFile);
end
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = raw{1};

mjdList  = zeros(0,1);
typeList = strings(0,1);
descList = strings(0,1);

stationLength = numel(station);
for k = 1:numel(lines)
    line = strtrim(lines{k});

    % Reject other stations before the relatively expensive STRSPLIT.
    if numel(line) <= stationLength || ...
            ~strcmpi(line(1:stationLength), station) || ...
            ~isspace(line(stationLength + 1))
        continue;
    end

    tok = strsplit(line);
    if numel(tok) < 4
        continue;
    end

    thisMjd = yymmmdd2mjd(tok{2});
    code = str2double(tok{3});

    if code == 1
        mjdList(end+1,1)  = thisMjd; %#ok<AGROW>
        typeList(end+1,1) = "EQUIP"; %#ok<AGROW>
        descList(end+1,1) = strrep(strjoin(tok(4:end), ' '), '_', ' '); %#ok<AGROW>
    elseif code == 2 && numel(tok) >= 7
        mag    = tok{6};
        eqid   = tok{7};
        mjdList(end+1,1)  = thisMjd; %#ok<AGROW>
        typeList(end+1,1) = "EQ"; %#ok<AGROW>
        descList(end+1,1) = sprintf('M%s %s', mag, eqid); %#ok<AGROW>
    end
end

MJD = mjdList;
Type = typeList;
Description = descList;
T = table(MJD, Type, Description);
if height(T) > 0
    T = sortrows(T, 'MJD');
end

end
