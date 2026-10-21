function A = offset_design(decimalYear, harmonicCount, offsetIndices)
% Builds the model matrix used to fit the time series.

t = decimalYear(:);
if isempty(t) || any(~isfinite(t))
    error('offset:InvalidEpochs', 'Decimal-year epochs must be finite and nonempty.');
end
if ~isscalar(harmonicCount) || ~isfinite(harmonicCount) || ...
        harmonicCount < 0 || harmonicCount ~= fix(harmonicCount)
    error('offset:InvalidHarmonicCount', ...
        'The harmonic count must be a nonnegative integer.');
end

t = t - t(1);
n = 2 + 2 .* harmonicCount + numel(offsetIndices);
A = zeros(numel(t), n);
A(:, 1:2) = [ones(numel(t), 1), t];

for k = 1:harmonicCount
    omegaT = 2 * pi * k * t;
    column = 2 + 2 .* (k - 1);
    A(:, column + (1:2)) = [cos(omegaT), sin(omegaT)];
end

for k = 1:numel(offsetIndices)
    idx = offsetIndices(k);
    if ~isscalar(idx) || ~isfinite(idx) || idx ~= fix(idx) || ...
            idx < 1 || idx > numel(t)
        error('offset:InvalidOffsetIndex', ...
            'Offset indices must be integer observation indices within the series.');
    end
    column = 2 + 2 .* harmonicCount + k;
    A(idx:end, column) = 1;
end

end
