function results = CommonMAD(time, data, k, momFile)
% COMMONMAD: One-pass MAD outlier detection on deterministic-model residuals.
%
% The function first fits a GNSS deterministic model:
%   bias + linear trend + annual + semiannual + offsets from MOM header
% Outliers are then detected once on the residuals using:
%   median(residuals) +/- k * 1.4826 * MAD(residuals)
%
% Parameters:
%   time    : timestamps, e.g. MJD
%   data    : coordinate values
%   k       : MAD tolerance coefficient, default 3
%   momFile : optional MOM file. Header lines "# offset MJD" are fitted.
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

if nargin < 3 || isempty(k) || isnan(k)
    k = 3;
end
if nargin < 4
    momFile = '';
end
if ~isvector(time) || ~isvector(data) || numel(time) ~= numel(data)
    error('TIME and DATA must be vectors of the same length.');
end

offsetTimes = readMomOffsets(momFile);
madScale = 1.4826;
originalShape = size(data);
t = time(:);
originalData = data(:);
cleanData = originalData;

if isempty(cleanData) || all(~isfinite(cleanData))
    warning('Input data is empty or contains no finite values.');
    results = makeResults(cleanData, originalData, originalShape, [], [], [], [], k, offsetTimes, momFile);
    return;
end

[residuals, fitted] = fitDeterministicModel(t, cleanData, offsetTimes);
validIdx = find(isfinite(t) & isfinite(cleanData) & isfinite(residuals));

if numel(validIdx) < 2
    warning('Not enough valid residuals (<2).');
    results = makeResults(cleanData, originalData, originalShape, [], [], residuals, fitted, k, offsetTimes, momFile);
    return;
end

residualValues = residuals(validIdx);
center = median(residualValues);
scale = madScale * median(abs(residualValues - center));

if scale == 0 || ~isfinite(scale)
    warning('Residual MAD is zero or invalid. No outliers detected.');
    results = makeResults(cleanData, originalData, originalShape, [], [], residuals, fitted, k, offsetTimes, momFile);
    return;
end

lowerBound = center - k * scale;
upperBound = center + k * scale;
isOutlier = residualValues < lowerBound | residualValues > upperBound;
removedIndices = validIdx(isOutlier);
removedTimes = t(removedIndices);
cleanData(removedIndices) = NaN;

results = makeResults(cleanData, originalData, originalShape, removedIndices, ...
    removedTimes, residuals, fitted, k, offsetTimes, momFile);
results.lowerBound = lowerBound;
results.upperBound = upperBound;
results.center = center;
results.scale = scale;
end

function results = makeResults(cleanData, originalData, originalShape, removedIndices, ...
    removedTimes, residuals, fitted, k, offsetTimes, momFile)
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
results.offsetTimes = offsetTimes(:);
results.momFile = momFile;
results.method = 'one-pass deterministic-model residual MAD';
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
