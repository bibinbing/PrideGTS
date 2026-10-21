function mask = offset_cand(m, accepted)
% Returns every epoch not already present in the offset model.

if ~isscalar(m) || ~isfinite(m) || m < 0 || m ~= fix(m)
    error('offset:InvalidObservationCount', ...
        'The observation count must be a nonnegative integer.');
end

mask = true(m, 1);
if isempty(accepted)
    return;
end

accepted = accepted(:);
if any(~isfinite(accepted)) || any(accepted ~= fix(accepted)) || ...
        any(accepted < 1) || any(accepted > m)
    error('offset:InvalidAcceptedOffsets', ...
        'Accepted offset indices must be integer observation indices.');
end
mask(unique(accepted)) = false;

end
