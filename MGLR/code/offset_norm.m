function opts = offset_norm(opts, projectRoot)
% Completes missing settings and prepares output folders.

if nargin < 2 || isempty(projectRoot)
    projectRoot = pwd;
end

defaults = offset_defaults();
opts = mergeStructs(defaults, opts);

validateScalar(opts.algorithm.alpha, 'algorithm.alpha');
if opts.algorithm.alpha <= 0 || opts.algorithm.alpha >= 1
    error('offset:InvalidAlpha', 'algorithm.alpha must be strictly between 0 and 1.');
end
validateInteger(opts.algorithm.harmonicCount, 'algorithm.harmonicCount', 0);
validateScalar(opts.algorithm.maxConditionNumber, 'algorithm.maxConditionNumber');
if opts.algorithm.maxConditionNumber <= 1
    error('offset:InvalidConditionLimit', ...
        'algorithm.maxConditionNumber must be greater than 1.');
end
thresholdMethod = char(lower(strtrim(string(opts.algorithm.thresholdMethod))));
if ~strcmp(thresholdMethod, 'station-bonferroni')
    error('offset:InvalidThresholdMethod', ...
        ['algorithm.thresholdMethod must be ''station-bonferroni'' for ', ...
        'the production offset workflow.']);
end
opts.algorithm.thresholdMethod = thresholdMethod;

validateScalar(opts.quality.transientMaxGapDays, ...
    'quality.transientMaxGapDays');
validateScalar(opts.quality.transientMaxDirectionCosine, ...
    'quality.transientMaxDirectionCosine');
validateScalar(opts.quality.transientMaxNetToGrossRatio, ...
    'quality.transientMaxNetToGrossRatio');
validateInteger(opts.quality.minSideObservations, ...
    'quality.minSideObservations', 1);
validateScalar(opts.quality.minSideSpanDays, ...
    'quality.minSideSpanDays');
if opts.quality.transientMaxGapDays <= 0 || ...
        opts.quality.transientMaxDirectionCosine < -1 || ...
        opts.quality.transientMaxDirectionCosine > 1 || ...
        opts.quality.transientMaxNetToGrossRatio < 0 || ...
        opts.quality.transientMaxNetToGrossRatio > 1 || ...
        opts.quality.minSideSpanDays <= 0
    error('offset:InvalidQualityControl', ...
        ['Transient and segment-support quality-control settings are ', ...
        'outside their valid ranges.']);
end

validateInteger(opts.noise.lsVceIterations, 'noise.lsVceIterations', 1);
validateScalar(opts.noise.lsVceTolerance, 'noise.lsVceTolerance');
if opts.noise.lsVceTolerance <= 0
    error('offset:InvalidLsVceTolerance', 'noise.lsVceTolerance must be positive.');
end
validateScalar(opts.noise.ridge, 'noise.ridge');
validateScalar(opts.noise.minVariance, 'noise.minVariance');
if opts.noise.ridge < 0 || opts.noise.minVariance <= 0
    error('offset:InvalidNoiseRegularization', ...
        'noise.ridge must be nonnegative and noise.minVariance must be positive.');
end

validateScalar(opts.figures.resolution, 'figures.resolution');
if opts.figures.resolution <= 0
    error('offset:InvalidFigureResolution', 'figures.resolution must be positive.');
end

if ischar(opts.inputPaths) || isstring(opts.inputPaths)
    opts.inputPaths = cellstr(opts.inputPaths);
end

if isempty(opts.outputDir)
    opts.outputDir = fullfile(projectRoot, 'result');
end
if ~isfolder(opts.outputDir)
    mkdir(opts.outputDir);
end

opts.output.seriesPath = fullfile(opts.outputDir, opts.output.seriesDir);
opts.output.figurePath = fullfile(opts.outputDir, opts.output.figureDir);

ensureDir(opts.output.seriesPath);
ensureDir(opts.output.figurePath);

end

function out = mergeStructs(base, override)
out = base;
if isempty(override)
    return;
end
names = fieldnames(override);
for i = 1:numel(names)
    name = names{i};
    if isstruct(override.(name)) && isfield(out, name) && isstruct(out.(name))
        out.(name) = mergeStructs(out.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end
end

function ensureDir(pathName)
if ~isfolder(pathName)
    mkdir(pathName);
end
end

function validateScalar(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
    error('offset:InvalidOption', '%s must be a finite real scalar.', name);
end
end

function validateInteger(value, name, lowerBound)
validateScalar(value, name);
if value ~= fix(value) || value < lowerBound
    error('offset:InvalidOption', '%s must be an integer greater than or equal to %d.', ...
        name, lowerBound);
end
end
