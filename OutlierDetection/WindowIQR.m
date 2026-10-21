function results = WindowIQR(time, data, windowSize, k, momFile)
% WINDOWIQR: One-pass windowed IQR outlier detection on model residuals.
%
% The function first fits a GNSS deterministic model:
%   bias + linear trend + annual + semiannual + offsets from MOM header
% Each epoch is then tested once against residuals in a centered moving
% window using:
%   median(window residuals) +/- k * IQR(window residuals)
%
% Parameters:
%   time       : timestamps, e.g. MJD
%   data       : coordinate values
%   windowSize : centered window length in samples, default 730
%   k          : IQR tolerance coefficient, default 3
%   momFile    : optional MOM file. Header lines "# offset MJD" are fitted.
%
% Returns:
%   results.data           : data with detected outliers replaced by NaN
%   results.removedIndices : indices of detected outliers
%   results.removedTimes   : timestamps of detected outliers
%   results.residuals      : deterministic-model residuals
%   results.fitted         : fitted deterministic model values
%   results.offsetTimes    : offset epochs read from the MOM header
%
% By Xianpeng Huang
% Updated to follow Hector residual-based outlier theory.

if nargin < 3 || isempty(windowSize) || isnan(windowSize)
    windowSize = 730;
end
if nargin < 4 || isempty(k) || isnan(k)
    k = 3;
end
if nargin < 5
    momFile = '';
end
if ~isvector(time) || ~isvector(data) || numel(time) ~= numel(data)
    error('TIME and DATA must be vectors of the same length.');
end

offsetTimes = readMomOffsets(momFile);
windowSize = max(1, floor(windowSize));
halfWindow = floor(windowSize / 2);
originalShape = size(data);
t = time(:);
originalData = data(:);
cleanData = originalData;
n = numel(cleanData);

if isempty(cleanData) || all(~isfinite(cleanData))
    warning('Input data is empty or contains no finite values.');
    results = makeResults(cleanData, originalData, originalShape, [], [], [], [], k, windowSize, offsetTimes, momFile);
    return;
end

[residuals, fitted] = fitDeterministicModel(t, cleanData, offsetTimes);
validGlobal = isfinite(t) & isfinite(cleanData) & isfinite(residuals);
removedMask = false(size(cleanData));

for i = 1:n
    if ~validGlobal(i)
        continue;
    end

    winStart = max(1, i - halfWindow);
    winEnd = min(n, i + halfWindow);
    winIdx = winStart:winEnd;
    winValid = winIdx(validGlobal(winIdx));

    if numel(winValid) < 4
        continue;
    end

    sortedResiduals = sort(residuals(winValid));
    q1 = hectorQuantile(sortedResiduals, 0.25);
    q3 = hectorQuantile(sortedResiduals, 0.75);
    center = hectorQuantile(sortedResiduals, 0.50);
    scale = q3 - q1;

    if scale == 0 || ~isfinite(scale)
        continue;
    end

    lowerBound = center - k * scale;
    upperBound = center + k * scale;
    if residuals(i) < lowerBound || residuals(i) > upperBound
        removedMask(i) = true;
    end
end

removedIndices = find(removedMask);
removedTimes = t(removedIndices);
cleanData(removedIndices) = NaN;

results = makeResults(cleanData, originalData, originalShape, removedIndices, ...
    removedTimes, residuals, fitted, k, windowSize, offsetTimes, momFile);
end

function results = makeResults(cleanData, originalData, originalShape, removedIndices, ...
    removedTimes, residuals, fitted, k, windowSize, offsetTimes, momFile)
if isempty(residuals)
    residuals = NaN(size(originalData));
end
if isempty(fitted)
    fitted = NaN(size(originalData));
end
results = struct();
results.data = reshape(cleanData, originalShape);
results.removedIndices = removedIndices(:);
results.removedTimes = removedTimes(:);
results.residuals = reshape(residuals, originalShape);
results.fitted = reshape(fitted, originalShape);
results.k = k;
results.windowSize = windowSize;
results.offsetTimes = offsetTimes(:);
results.momFile = momFile;
results.method = 'one-pass deterministic-model residual windowed IQR';
results.model = 'bias + linear trend + annual + semiannual + MOM-header offsets';
end

function [residuals, fitted] = fitDeterministicModel(t, y, offsetTimes)
residuals = NaN(size(y));
fitted = NaN(size(y));
valid = isfinite(t) & isfinite(y);
if ~any(valid)
    return;
end

tv = t(valid);
yv = y(valid);
referenceEpoch = 0.5 * (min(tv) + max(tv));
H = buildDesignMatrix(tv, referenceEpoch, offsetTimes);

theta = pinv(H) * yv;
finiteTime = isfinite(t);
Hfull = buildDesignMatrix(t(finiteTime), referenceEpoch, offsetTimes);
fitted(finiteTime) = Hfull * theta;
residuals(valid) = y(valid) - fitted(valid);
end

function H = buildDesignMatrix(t, referenceEpoch, offsetTimes)
t = t(:);
phaseAnnual = 2 * pi * (t - 51544.0) / 365.25;
phaseSemiannual = 2 * pi * (t - 51544.0) / 182.625;
H = [ones(numel(t), 1), ...
     t - referenceEpoch, ...
     cos(phaseAnnual), sin(phaseAnnual), ...
     cos(phaseSemiannual), sin(phaseSemiannual)];

for i = 1:numel(offsetTimes)
    H(:, end + 1) = double(t + 1.0e-7 > offsetTimes(i)); %#ok<AGROW>
end
end

function offsetTimes = readMomOffsets(momFile)
offsetTimes = [];
if isempty(momFile)
    return;
end
if isstring(momFile)
    momFile = char(momFile);
end
if ~ischar(momFile)
    error('MOMFILE must be a character vector or string scalar.');
end

fid = fopen(momFile, 'r');
if fid < 0
    error('Could not open MOM file: %s', momFile);
end
cleaner = onCleanup(@() fclose(fid));

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    lineTrim = strtrim(line);
    if isempty(lineTrim)
        continue;
    end
    if lineTrim(1) ~= '#'
        break;
    end

    token = regexpi(lineTrim, '^#\s*offset\s+(.+)$', 'tokens', 'once');
    if ~isempty(token)
        value = sscanf(token{1}, '%f', 1);
        if ~isempty(value) && isfinite(value)
            offsetTimes(end + 1, 1) = value; %#ok<AGROW>
        end
    end
end
end

function value = hectorQuantile(sortedValues, fraction)
n = numel(sortedValues);
idx = floor(fraction * n) + 1;
idx = max(1, min(n, idx));
value = sortedValues(idx);
end
