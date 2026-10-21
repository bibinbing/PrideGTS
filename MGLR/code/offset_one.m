function result = offset_one(inputFile, opts)
% Runs the full offset workflow for one input file.

ts = offset_read(inputFile);
detection = offset_detect(ts, opts);
correctedNeu = offset_correct(ts.neu, detection);

outFile = offset_write(ts, correctedNeu, opts, detection);

figureFiles = {};
if opts.figures.enable
    figureFiles = offset_figs(ts, detection, correctedNeu, opts);
end

result = struct();
result.station = ts.station;
result.format = ts.sourceFormat;
result.outputFile = outFile;
result.figureFiles = figureFiles;
result.offsetTable = detection.offsetTable;
result.detection = detection;
result.timeSeries = ts;
result.correctedNeu = correctedNeu;
if isfield(ts, 'quality')
    result.quality = ts.quality;
else
    result.quality = struct();
end

end
