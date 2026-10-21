function [Qf, eigenvalues, eigenvectors] = offset_flicker(epochDays, usePersistentCache)
% Builds flicker-noise covariance while preserving gaps in the daily grid.

persistent cachedDayIndex cachedQf cachedEigenvalues cachedEigenvectors
if nargin < 2
    usePersistentCache = true;
end
if (ischar(epochDays) || (isstring(epochDays) && isscalar(epochDays))) && ...
        strcmpi(strtrim(string(epochDays)), "clear")
    cachedDayIndex = [];
    cachedQf = [];
    cachedEigenvalues = [];
    cachedEigenvectors = [];
    Qf = [];
    eigenvalues = [];
    eigenvectors = [];
    return;
end

dayIndex = normalizeDayIndex(epochDays);
m = numel(dayIndex);
if m == 0
    Qf = zeros(0, 0);
    eigenvalues = zeros(0, 1);
    eigenvectors = zeros(0, 0);
    return;
end

if usePersistentCache && ~isempty(cachedDayIndex) && ...
        isequal(dayIndex, cachedDayIndex)
    Qf = cachedQf;
    if nargout > 1
        if isempty(cachedEigenvalues) || ...
                (nargout > 2 && isempty(cachedEigenvectors))
            if nargout > 2
                [cachedEigenvectors, cachedEigenvalues] = eig(Qf, 'vector');
            else
                cachedEigenvalues = eig(Qf, 'vector');
            end
            validateEigenvalues(cachedEigenvalues);
        end
        eigenvalues = cachedEigenvalues;
    end
    if nargout > 2
        eigenvectors = cachedEigenvectors;
    end
    return;
end

dayCount = dayIndex(end) + 1;
d = 0.5;
h = zeros(dayCount, 1);
h(1) = 1;
for n = 1:dayCount - 1
    h(n + 1) = h(n) .* (n - 1 + d) ./ n;
end

% Maps full-grid days back to observed rows with a compact index.
observedRow = zeros(dayCount, 1, 'uint32');
observedRow(dayIndex + 1) = uint32(1:m);
Qf = zeros(m, m);

for lag = 0:dayCount - 1
    n = dayCount - lag;
    covarianceByEarlierDay = cumsum(h(1:n) .* h(1 + lag:dayCount));
    earlierRow = observedRow(1:n);
    laterRow = observedRow(1 + lag:dayCount);
    keep = earlierRow ~= 0 & laterRow ~= 0;
    if ~any(keep)
        continue;
    end

    cols = double(earlierRow(keep));
    rows = double(laterRow(keep));
    values = covarianceByEarlierDay(keep);
    Qf(sub2ind([m, m], rows, cols)) = values;
    if lag ~= 0
        Qf(sub2ind([m, m], cols, rows)) = values;
    end
end

if nargout > 1
    if nargout > 2
        [eigenvectors, eigenvalues] = eig(Qf, 'vector');
    else
        eigenvalues = eig(Qf, 'vector');
    end
    validateEigenvalues(eigenvalues);
end
if usePersistentCache
    cachedDayIndex = dayIndex;
    cachedQf = Qf;
    if nargout > 1
        cachedEigenvalues = eigenvalues;
    else
        cachedEigenvalues = [];
    end
    if nargout > 2
        cachedEigenvectors = eigenvectors;
    else
        cachedEigenvectors = [];
    end
end

end

function validateEigenvalues(eigenvalues)
if any(~isfinite(eigenvalues)) || any(~isreal(eigenvalues)) || ...
        any(eigenvalues <= 0)
    error('offset:InvalidFlickerSpectrum', ...
        'The Hosking covariance eigenspectrum must be finite and positive.');
end
end

function dayIndex = normalizeDayIndex(epochDays)
if isempty(epochDays)
    dayIndex = zeros(0, 1);
    return;
end

if isscalar(epochDays)
    m = double(epochDays);
    if ~isfinite(m) || m < 0 || m ~= fix(m)
        error('offset:InvalidFlickerSize', ...
            'A scalar flicker-noise size must be a nonnegative integer.');
    end
    dayIndex = (0:m - 1).';
    return;
end

epochs = double(epochDays(:));
if any(~isfinite(epochs)) || any(diff(epochs) <= 0)
    error('offset:InvalidFlickerEpochs', ...
        'Flicker-noise epochs must be finite, unique, and strictly increasing.');
end

relativeDays = epochs - epochs(1);
dayIndex = round(relativeDays);
roundingTolerance = 1e-7 .* max(1, max(abs(relativeDays)));
if any(abs(relativeDays - dayIndex) > roundingTolerance)
    error('offset:NonDailyFlickerEpochs', ...
        ['The Hosking covariance requires epochs on a daily grid. ', ...
        'MJD values may contain missing days but not fractional-day spacing.']);
end
if any(diff(dayIndex) <= 0)
    error('offset:DuplicateFlickerDays', ...
        'Epochs map to duplicate days in the Hosking covariance.');
end

end
