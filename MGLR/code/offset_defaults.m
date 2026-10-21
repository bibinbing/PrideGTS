function opts = offset_defaults()
% Provides default settings for the offset workflow.

opts = struct();
opts.inputPaths = {};
opts.outputDir = fullfile(pwd, 'result');

opts.algorithm = struct();
opts.algorithm.alpha = 0.001;
opts.algorithm.harmonicCount = 2;
opts.algorithm.maxConditionNumber = 1e12;

% Uses a station-level threshold across the current candidates.
opts.algorithm.thresholdMethod = 'station-bonferroni';

opts.quality = struct();
opts.quality.transientMaxGapDays = 7;
opts.quality.transientMaxDirectionCosine = -0.8;
opts.quality.transientMaxNetToGrossRatio = 0.2;
opts.quality.minSideObservations = 30;
opts.quality.minSideSpanDays = 30;

opts.noise = struct();
opts.noise.model = 'white-flicker';
opts.noise.lsVceIterations = 20;
opts.noise.lsVceTolerance = 1e-3;
opts.noise.ridge = 1e-10;
opts.noise.minVariance = 1e-12;

opts.figures = struct();
opts.figures.enable = true;
opts.figures.resolution = 300;

opts.output = struct();
opts.output.seriesDir = 'timeSeries';
opts.output.figureDir = 'fig';

end
