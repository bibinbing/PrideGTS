function figureFiles = offset_figs(ts, detection, correctedNeu, opts)
% Saves the standard output figures for one station.

stem = offset_safe(ts.station);

comparisonFile = fullfile(opts.output.figurePath, [stem '_comparison.png']);
powerFile = fullfile(opts.output.figurePath, [stem '_offset_power.png']);
residualFile = fullfile(opts.output.figurePath, [stem '_residuals.png']);
vectorFile = fullfile(opts.output.figurePath, [stem '_offset_vectors.png']);

comparisonOut = offset_plotCompare(ts, detection, correctedNeu, comparisonFile, opts);
powerOut = offset_plotPower(ts, detection, powerFile, opts);
residualOut = offset_plotResiduals(ts, detection, residualFile, opts);
vectorOut = offset_plotVectors(ts, detection, vectorFile, opts);

figureFiles = {comparisonOut; powerOut; residualOut; vectorOut};

end
