function threshold = offset_stationThreshold(alpha, componentCount, candidateCount)
% Returns the Bonferroni station-level threshold for the current GLR scan.

if ~isscalar(candidateCount) || ~isfinite(candidateCount) || ...
        candidateCount < 1 || candidateCount ~= fix(candidateCount)
    error('offset:InvalidCandidateCount', ...
        'The identifiable candidate count must be a positive integer.');
end

threshold = offset_chi2crit(alpha ./ candidateCount, componentCount);

end
