function correctedNeu = offset_correct(neu, detection)
% Removes detected offsets from the NEU time series.

correctedNeu = neu;
if isfield(detection, 'appliedIndices') && ...
        isfield(detection, 'appliedCoefficients')
    indices = detection.appliedIndices(:);
    coefficients = detection.appliedCoefficients;
else
    indices = detection.acceptedIndices(:);
    coefficients = detection.offsetCoefficients;
end
if ~isequal(size(coefficients), [numel(indices), size(neu, 2)])
    error('offset:CorrectionAlignment', ...
        'Correction epochs and component estimates are not aligned.');
end
for i = 1:numel(indices)
    idx = indices(i);
    correctedNeu(idx:end, :) = correctedNeu(idx:end, :) - coefficients(i, :);
end

end
