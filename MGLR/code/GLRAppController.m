classdef GLRAppController < handle
    % Coordinates station import, offset detection, logging, and plotting.

    properties (Access = private)
        App
        ProjectRoot
        Stations = struct([])
        Details = struct([])
        CurrentDetailIndex = []
        LogLines = cell(0, 1)
        IsBusy = false
        ProgressDialog = []
    end

    methods
        function obj = GLRAppController(app, projectRoot)
            obj.App = app;
            obj.ProjectRoot = projectRoot;
            addpath(projectRoot);
            addpath(fullfile(projectRoot, 'code'));
            obj.initializeUi();
        end

        function browseInput(obj)
            obj.resetStateButton(obj.App.DataImportBrowseButton);
            startFolder = obj.existingFolder( ...
                obj.App.DataImportEditField.Value, ...
                fullfile(obj.ProjectRoot, 'data'));
            selectedFolder = uigetdir(startFolder, ...
                'Select a folder containing station MOM files');
            if isequal(selectedFolder, 0)
                return;
            end
            obj.App.DataImportEditField.Value = selectedFolder;
            obj.refreshStationTree(true);
        end

        function browseOutput(obj)
            obj.resetStateButton(obj.App.DataOutputBrowseButton);
            startFolder = obj.existingFolder( ...
                obj.App.DataOutputEditField.Value, ...
                fullfile(obj.ProjectRoot, 'result'));
            selectedFolder = uigetdir(startFolder, 'Select output folder');
            if ~isequal(selectedFolder, 0)
                obj.App.DataOutputEditField.Value = selectedFolder;
            end
        end

        function refreshStationTree(obj, writeLog)
            if nargin < 2
                writeLog = true;
            end
            try
                inputText = strtrim(string(obj.App.DataImportEditField.Value));
                if strlength(inputText) == 0
                    inputText = string(fullfile(obj.ProjectRoot, 'data'));
                    obj.App.DataImportEditField.Value = char(inputText);
                end
                inputFiles = offset_inputs(obj.parseInputPaths(inputText));
                obj.rebuildStationTree(inputFiles);
                if writeLog
                    obj.appendLog(sprintf( ...
                        'Imported %d station group(s) from %s.', ...
                        numel(inputFiles), char(inputText)));
                end
            catch ME
                obj.Stations = struct([]);
                delete(obj.App.Tree.Children);
                obj.appendLog(['Station import failed: ' ME.message]);
                uialert(obj.App.GLRUIFigure, ME.message, 'Station Import');
            end
        end

        function runAll(obj)
            obj.resetStateButton(obj.App.RunallButton);
            if obj.IsBusy
                return;
            end
            obj.refreshStationTree(false);
            if isempty(obj.Stations)
                uialert(obj.App.GLRUIFigure, ...
                    'No complete station MOM groups were found.', 'Run All');
                return;
            end
            inputFiles = {obj.Stations.inputFile};
            obj.runInputs(inputFiles, 'Run all');
        end

        function runSelected(obj)
            obj.resetStateButton(obj.App.RunSelectedButton);
            if obj.IsBusy
                return;
            end
            selectedNodes = obj.App.Tree.SelectedNodes;
            if isempty(selectedNodes)
                uialert(obj.App.GLRUIFigure, ...
                    'Select a station or one of its MOM component files first.', ...
                    'Run Selected');
                return;
            end
            nodeData = selectedNodes(1).NodeData;
            if ~isstruct(nodeData) || ~isfield(nodeData, 'inputFile')
                uialert(obj.App.GLRUIFigure, ...
                    'The selected tree item is not a station MOM group.', ...
                    'Run Selected');
                return;
            end
            obj.runInputs({nodeData.inputFile}, ...
                sprintf('Run selected: %s', nodeData.station));
        end

        function exportCurrentFigure(obj)
            obj.resetStateButton(obj.App.ExportFigureButton);
            if obj.IsBusy
                return;
            end
            if isempty(obj.CurrentDetailIndex) || isempty(obj.Details)
                uialert(obj.App.GLRUIFigure, ...
                    'Run offset detection and select a successful station first.', ...
                    'Export Figure');
                return;
            end
            detail = obj.Details(obj.CurrentDetailIndex);
            if ~detail.success || ~isfile(detail.cacheFile)
                uialert(obj.App.GLRUIFigure, ...
                    'The selected station has no result figures to export.', ...
                    'Export Figure');
                return;
            end

            outputFolder = uigetdir(obj.existingFolder( ...
                obj.App.DataOutputEditField.Value, ...
                fullfile(obj.ProjectRoot, 'result')), ...
                'Select a folder for the four result figures');
            if isequal(outputFolder, 0)
                return;
            end
            try
                plotDetail = obj.loadPlotDetail(detail.cacheFile);
                opts = offset_defaults();
                stem = offset_safe(plotDetail.station);
                outputFiles = {
                    fullfile(outputFolder, [stem '_time_series.png'])
                    fullfile(outputFolder, [stem '_residuals.png'])
                    fullfile(outputFolder, [stem '_glr_power.png'])
                    fullfile(outputFolder, [stem '_step_vector.png'])};
                offset_plotCompare(plotDetail.ts, plotDetail.detection, ...
                    plotDetail.correctedNeu, outputFiles{1}, opts);
                offset_plotResiduals(plotDetail.ts, plotDetail.detection, ...
                    outputFiles{2}, opts);
                offset_plotPower(plotDetail.ts, plotDetail.detection, ...
                    outputFiles{3}, opts);
                offset_plotVectors(plotDetail.ts, plotDetail.detection, ...
                    outputFiles{4}, opts);
                obj.appendLog(sprintf( ...
                    '%s: exported four result figures to %s.', ...
                    plotDetail.station, outputFolder));
                uialert(obj.App.GLRUIFigure, sprintf( ...
                    'Four figures for %s were saved to:\n%s', ...
                    plotDetail.station, outputFolder), ...
                    'Export Figure', 'Icon', 'success');
            catch ME
                obj.appendLog(['Figure export failed: ' ME.message]);
                uialert(obj.App.GLRUIFigure, ME.message, 'Export Figure');
            end
        end
    end

    methods (Access = private)
        function initializeUi(obj)
            if strlength(strtrim(string( ...
                    obj.App.DataImportEditField.Value))) == 0
                obj.App.DataImportEditField.Value = ...
                    fullfile(obj.ProjectRoot, 'data');
            end
            if strlength(strtrim(string( ...
                    obj.App.DataOutputEditField.Value))) == 0
                obj.App.DataOutputEditField.Value = ...
                    fullfile(obj.ProjectRoot, 'result');
            end
            obj.App.TextArea.Editable = 'off';
            obj.App.TextArea.Value = {''};
            obj.App.Tree.Multiselect = 'off';
            obj.App.Tree.SelectionChangedFcn = ...
                @(~, event) obj.treeSelectionChanged(event);
            obj.App.DataImportEditField.ValueChangedFcn = ...
                @(~, ~) obj.refreshStationTree(true);
            obj.appendLog(['Ready. Import MOM files, select a station, ', ...
                'then run joint 3-D GLR offset detection.']);
            obj.refreshStationTree(false);
            obj.clearPlots();
        end

        function paths = parseInputPaths(~, rawValue)
            pieces = strtrim(split(string(rawValue), ';'));
            pieces(pieces == "") = [];
            paths = cellstr(pieces);
        end

        function folder = existingFolder(~, rawValue, fallback)
            folder = fallback;
            if ~isfolder(folder)
                parent = fileparts(folder);
                if isfolder(parent)
                    folder = parent;
                else
                    folder = pwd;
                end
            end
            candidates = strtrim(split(string(rawValue), ';'));
            for i = 1:numel(candidates)
                candidate = char(candidates(i));
                if isfolder(candidate)
                    folder = candidate;
                    return;
                end
                parent = fileparts(candidate);
                if isfolder(parent)
                    folder = parent;
                    return;
                end
            end
        end

        function rebuildStationTree(obj, inputFiles)
            delete(obj.App.Tree.Children);
            obj.Stations = repmat(obj.emptyStation(), numel(inputFiles), 1);
            componentNames = {'North', 'East', 'Up'};
            for i = 1:numel(inputFiles)
                inputFile = inputFiles{i};
                [folder, stem, extension] = fileparts(inputFile);
                token = regexp(stem, '^(.*)_[012]$', 'tokens', 'once');
                if isempty(token)
                    base = stem;
                else
                    base = token{1};
                end
                station = upper(char(base));
                componentFiles = cell(1, 3);
                nodeData = struct('station', station, ...
                    'inputFile', inputFile);
                stationNode = uitreenode(obj.App.Tree, ...
                    'Text', station, 'NodeData', nodeData);
                for componentIndex = 1:3
                    componentFiles{componentIndex} = fullfile(folder, ...
                        sprintf('%s_%d%s', base, componentIndex - 1, extension));
                    [~, fileName, fileExtension] = ...
                        fileparts(componentFiles{componentIndex});
                    label = sprintf('%s - %s%s', ...
                        componentNames{componentIndex}, fileName, fileExtension);
                    if ~isfile(componentFiles{componentIndex})
                        label = [label ' [missing]'];
                    end
                    uitreenode(stationNode, 'Text', label, ...
                        'NodeData', nodeData);
                end
                obj.Stations(i).station = station;
                obj.Stations(i).inputFile = inputFile;
                obj.Stations(i).componentFiles = componentFiles;
                obj.Stations(i).node = stationNode;
            end
            if ~isempty(obj.App.Tree.Children)
                expand(obj.App.Tree, 'all');
            end
        end

        function station = emptyStation(~)
            station = struct('station', '', 'inputFile', '', ...
                'componentFiles', {{}}, 'node', []);
        end

        function detail = emptyDetail(~)
            detail = struct('success', false, 'station', '', ...
                'inputFile', '', 'outputFile', '', 'cacheFile', '', ...
                'message', '', 'counts', struct());
        end

        function counts = emptyCounts(~)
            counts = struct('total', 0, 'known', 0, 'persistent', 0, ...
                'lowSupport', 0, 'transientEvents', 0, ...
                'transientPairs', 0, 'applied', 0, ...
                'glrCandidates', 0, 'failed', 0);
        end

        function runInputs(obj, inputFiles, runName)
            obj.IsBusy = true;
            previousTitle = obj.App.GLRUIFigure.Name;
            obj.setControlsBusy(true);
            cleanup = onCleanup(@() obj.restoreUi(previousTitle));
            total = numel(inputFiles);
            aggregate = obj.emptyCounts();
            successCount = 0;

            try
                opts = obj.buildOptions();
                plotDataFolder = fullfile(opts.outputDir, 'plotData');
                if ~isfolder(plotDataFolder)
                    mkdir(plotDataFolder);
                end
                obj.Details = repmat(obj.emptyDetail(), total, 1);
                obj.CurrentDetailIndex = [];
                obj.clearPlots();
                obj.appendLog({sprintf('%s started for %d station(s).', ...
                    runName, total), ...
                    ['Offset types: known-constrained = MOM header constraint; ', ...
                    'persistent-candidate = supported GLR step; ', ...
                    'low-support-candidate = insufficient side support; ', ...
                    'transient-pair = close cancelling pair.']});

                obj.ProgressDialog = uiprogressdlg(obj.App.GLRUIFigure, ...
                    'Title', 'GLR offset detection', ...
                    'Message', sprintf('Preparing %d station(s)...', total), ...
                    'Value', 0);
                obj.App.GLRUIFigure.Name = sprintf( ...
                    'GLR - processing %d station(s)', total);

                for i = 1:total
                    inputFile = inputFiles{i};
                    station = obj.stationName(inputFile);
                    obj.updateStationNode(inputFile, 'running');
                    obj.ProgressDialog.Value = (i - 1) / total;
                    obj.ProgressDialog.Message = sprintf( ...
                        'Processing %d of %d: %s', i, total, station);
                    obj.appendLog(sprintf('%s: detection started.', station));
                    drawnow;
                    stationTimer = tic;

                    detail = obj.emptyDetail();
                    detail.station = station;
                    detail.inputFile = inputFile;
                    try
                        result = obj.runOneResponsively(inputFile, opts);
                        [summary, eventLines, warningLines, counts] = ...
                            obj.summarizeResult(result);
                        detail.success = true;
                        detail.station = result.station;
                        detail.outputFile = result.outputFile;
                        detail.cacheFile = obj.writePlotCache( ...
                            result, plotDataFolder, i);
                        detail.message = summary;
                        detail.counts = counts;
                        obj.Details(i) = detail;
                        successCount = successCount + 1;
                        aggregate = obj.addCounts(aggregate, counts);
                        obj.appendLog(sprintf('%s: %s; elapsed %.1f s.', ...
                            result.station, summary, toc(stationTimer)));
                        obj.appendLog(eventLines);
                        obj.appendLog(warningLines);
                        obj.updateStationNode(inputFile, sprintf( ...
                            'completed: %d offset(s)', counts.total));
                        clear result
                    catch stationError
                        detail.message = stationError.message;
                        obj.Details(i) = detail;
                        aggregate.failed = aggregate.failed + 1;
                        obj.updateStationNode(inputFile, 'failed');
                        obj.appendLog(sprintf( ...
                            '%s: failed after %.1f s - %s', ...
                            station, toc(stationTimer), stationError.message));
                    end
                    obj.ProgressDialog.Value = i / total;
                    drawnow;
                end

                summary = obj.aggregateSummary( ...
                    runName, successCount, total, aggregate);
                obj.appendLog(summary);
                firstSuccess = find([obj.Details.success], 1, 'first');
                if isempty(firstSuccess)
                    obj.App.GLRUIFigure.Name = 'GLR - analysis failed';
                    uialert(obj.App.GLRUIFigure, ...
                        'No selected station completed successfully. See Status & Log.', ...
                        'GLR Offset Detection');
                else
                    obj.showDetail(firstSuccess);
                    obj.App.GLRUIFigure.Name = sprintf( ...
                        'GLR - completed %d/%d', successCount, total);
                    obj.ProgressDialog.Message = summary;
                end
            catch ME
                obj.App.GLRUIFigure.Name = 'GLR - run failed';
                obj.appendLog(['Run failed: ' ME.message]);
                uialert(obj.App.GLRUIFigure, ME.message, 'GLR Offset Detection');
            end
            clear cleanup
        end

        function opts = buildOptions(obj)
            outputText = strtrim(string(obj.App.DataOutputEditField.Value));
            if strlength(outputText) == 0
                outputText = string(fullfile(obj.ProjectRoot, 'result'));
                obj.App.DataOutputEditField.Value = char(outputText);
            end
            opts = offset_defaults();
            opts.outputDir = char(outputText);
            opts.figures.enable = false;
            opts = offset_norm(opts, obj.ProjectRoot);
        end

        function result = runOneResponsively(obj, inputFile, opts)
            future = [];
            try
                workerOpts = opts;
                if isfield(workerOpts, 'runtime') && ...
                        isstruct(workerOpts.runtime)
                    workerOpts.runtime.manageThreads = false;
                end
                future = parfeval(backgroundPool, @offset_one, 1, ...
                    inputFile, workerOpts);
            catch
                % Foreground execution is used if backgroundPool is unavailable.
            end
            if isempty(future)
                result = offset_one(inputFile, opts);
                return;
            end
            futureCleanup = onCleanup(@() obj.cancelFuture(future));
            while ~strcmp(future.State, 'finished')
                if isempty(obj.App.GLRUIFigure) || ...
                        ~isvalid(obj.App.GLRUIFigure)
                    error('offset:AppClosed', ...
                        'The app was closed before detection finished.');
                end
                drawnow limitrate;
                pause(0.05);
            end
            if ~isempty(future.Error)
                error('offset:BackgroundFailure', '%s', future.Error.message);
            end
            result = fetchOutputs(future);
            clear futureCleanup
        end

        function cancelFuture(~, future)
            try
                if ~isempty(future) && isvalid(future) && ...
                        ~strcmp(future.State, 'finished')
                    cancel(future);
                end
            catch
                % Cancellation is best-effort during cleanup.
            end
        end

        function cacheFile = writePlotCache(~, result, folder, row)
            plotDetail = struct('station', result.station, ...
                'ts', result.timeSeries, 'detection', result.detection, ...
                'correctedNeu', result.correctedNeu);
            cacheFile = fullfile(folder, sprintf('%03d_%s_plot_data.mat', ...
                row, offset_safe(result.station)));
            save(cacheFile, 'plotDetail', '-v7.3');
        end

        function plotDetail = loadPlotDetail(~, cacheFile)
            if ~isfile(cacheFile)
                error('offset:MissingPlotData', ...
                    'Plot data file was not found: %s', cacheFile);
            end
            payload = load(cacheFile, 'plotDetail');
            if ~isfield(payload, 'plotDetail')
                error('offset:InvalidPlotData', ...
                    'Plot data are invalid: %s', cacheFile);
            end
            plotDetail = payload.plotDetail;
        end

        function [summary, eventLines, warningLines, counts] = ...
                summarizeResult(obj, result)
            tableData = result.offsetTable;
            counts = obj.emptyCounts();
            counts.total = height(tableData);
            if counts.total > 0
                classes = string(tableData.EventClass);
                counts.known = nnz(classes == "known-constrained");
                counts.persistent = nnz(classes == "persistent-candidate");
                counts.lowSupport = nnz(classes == "low-support-candidate");
                counts.transientEvents = nnz(classes == "transient-pair");
                counts.transientPairs = numel(unique( ...
                    tableData.PairID(tableData.PairID > 0)));
                counts.applied = nnz(tableData.ApplyCorrection);
            end
            if isfield(result.detection, 'candidateCount')
                counts.glrCandidates = result.detection.candidateCount;
            end
            summary = sprintf([ ...
                'total %d offset event(s): known-constrained %d, ', ...
                'persistent-candidate %d, low-support-candidate %d, ', ...
                'transient-pair %d event(s) / %d pair(s); ', ...
                'applied corrections %d; GLR candidates %d'], ...
                counts.total, counts.known, counts.persistent, ...
                counts.lowSupport, counts.transientEvents, ...
                counts.transientPairs, counts.applied, ...
                counts.glrCandidates);

            eventLines = cell(0, 1);
            for i = 1:height(tableData)
                if isfinite(tableData.Power(i))
                    powerText = sprintf('GLR %.3f / threshold %.3f', ...
                        tableData.Power(i), tableData.Threshold(i));
                else
                    powerText = 'GLR n/a (known constraint)';
                end
                if tableData.ApplyCorrection(i)
                    appliedText = 'applied';
                else
                    appliedText = 'not applied';
                end
                eventLines{end + 1, 1} = sprintf([ ...
                    '%s offset %d/%d | %s | %s | ', ...
                    'N %.2f, E %.2f, U %.2f mm | magnitude %.2f mm | ', ...
                    '%s | %s | quality %s'], ...
                    result.station, i, height(tableData), ...
                    char(string(tableData.Date(i))), ...
                    char(string(tableData.EventClass(i))), ...
                    tableData.North_mm(i), tableData.East_mm(i), ...
                    tableData.Up_mm(i), tableData.Magnitude_mm(i), ...
                    powerText, appliedText, ...
                    char(string(tableData.QualityFlag(i)))); %#ok<AGROW>
            end

            warningLines = cell(0, 1);
            if isfield(result.detection, 'warnings') && ...
                    ~isempty(result.detection.warnings)
                rawWarnings = result.detection.warnings;
                if ischar(rawWarnings)
                    rawWarnings = {rawWarnings};
                elseif isstring(rawWarnings)
                    rawWarnings = cellstr(rawWarnings(:));
                end
                if iscell(rawWarnings)
                    for i = 1:numel(rawWarnings)
                        warningLines{end + 1, 1} = sprintf( ...
                            '%s warning: %s', result.station, ...
                            char(string(rawWarnings{i}))); %#ok<AGROW>
                    end
                end
            end
            if isfield(result, 'quality') && isstruct(result.quality)
                removed = obj.scalarField( ...
                    result.quality, 'removedNonfiniteCount', 0);
                gaps = obj.scalarField(result.quality, 'gapCount', 0);
                warningLines{end + 1, 1} = sprintf( ...
                    '%s data quality: removed %d non-finite row(s); %d sampling gap(s).', ...
                    result.station, round(removed), round(gaps));
            end
            if isfield(result.detection, 'finalNoise')
                noise = result.detection.finalNoise;
                converged = obj.scalarField(noise, 'lsVceConverged', NaN);
                iterations = obj.scalarField(noise, 'lsVceIterations', NaN);
                warningLines{end + 1, 1} = sprintf( ...
                    '%s noise model: %s; LS-VCE converged %s; iterations %s.', ...
                    result.station, obj.textField(noise, 'model', 'unknown'), ...
                    obj.logicalText(converged), obj.numberText(iterations));
            end
        end

        function value = scalarField(~, source, fieldName, fallback)
            value = fallback;
            if isstruct(source) && isfield(source, fieldName)
                candidate = source.(fieldName);
                if isnumeric(candidate) || islogical(candidate)
                    if isscalar(candidate) && isfinite(double(candidate))
                        value = double(candidate);
                    end
                end
            end
        end

        function value = textField(~, source, fieldName, fallback)
            value = fallback;
            if isstruct(source) && isfield(source, fieldName) && ...
                    ~isempty(source.(fieldName))
                value = char(string(source.(fieldName)));
            end
        end

        function textValue = logicalText(~, value)
            if ~isfinite(value)
                textValue = 'n/a';
            elseif logical(value)
                textValue = 'yes';
            else
                textValue = 'no';
            end
        end

        function textValue = numberText(~, value)
            if isfinite(value)
                textValue = sprintf('%d', round(value));
            else
                textValue = 'n/a';
            end
        end

        function total = addCounts(~, total, value)
            fields = fieldnames(total);
            for i = 1:numel(fields)
                name = fields{i};
                if isfield(value, name)
                    total.(name) = total.(name) + value.(name);
                end
            end
        end

        function summary = aggregateSummary(~, runName, succeeded, total, counts)
            summary = sprintf([ ...
                '%s completed: %d/%d station(s) succeeded; ', ...
                'total offset events %d; known-constrained %d; ', ...
                'persistent-candidate %d; low-support-candidate %d; ', ...
                'transient-pair %d event(s) / %d pair(s); ', ...
                'applied corrections %d; GLR candidates %d; failures %d.'], ...
                runName, succeeded, total, counts.total, counts.known, ...
                counts.persistent, counts.lowSupport, ...
                counts.transientEvents, counts.transientPairs, ...
                counts.applied, counts.glrCandidates, counts.failed);
        end

        function treeSelectionChanged(obj, event)
            if isempty(event.SelectedNodes)
                return;
            end
            data = event.SelectedNodes(1).NodeData;
            if ~isstruct(data) || ~isfield(data, 'inputFile')
                return;
            end
            if isempty(obj.Details)
                obj.App.GLRUIFigure.Name = sprintf('GLR - %s selected', data.station);
                return;
            end
            index = find(strcmpi({obj.Details.inputFile}, data.inputFile), ...
                1, 'first');
            if isempty(index)
                index = find(strcmpi({obj.Details.station}, data.station), ...
                    1, 'first');
            end
            if ~isempty(index) && obj.Details(index).success
                obj.showDetail(index);
            end
        end

        function showDetail(obj, index)
            if index < 1 || index > numel(obj.Details) || ...
                    ~obj.Details(index).success
                return;
            end
            try
                plotDetail = obj.loadPlotDetail(obj.Details(index).cacheFile);
                obj.plotAll(plotDetail);
                obj.CurrentDetailIndex = index;
                obj.App.GLRUIFigure.Name = sprintf( ...
                    'GLR - %s', plotDetail.station);
            catch ME
                obj.appendLog(['Result display failed: ' ME.message]);
                uialert(obj.App.GLRUIFigure, ME.message, 'Display Result');
            end
        end

        function plotAll(obj, detail)
            obj.plotTimeSeries(detail);
            obj.plotResidual(detail);
            obj.plotGlrPower(detail);
            obj.plotStepVector(detail);
            drawnow limitrate;
        end

        function plotTimeSeries(obj, detail)
            axesList = [obj.App.NorthUIAxes, obj.App.EastUIAxes, ...
                obj.App.UpUIAxes];
            names = {'North', 'East', 'Up'};
            shortNames = {'N', 'E', 'U'};
            x = detail.ts.decimalYear(:);
            raw = 1000 .* detail.ts.neu;
            corrected = 1000 .* detail.correctedNeu;
            for component = 1:3
                ax = axesList(component);
                obj.resetAxes(ax);
                hold(ax, 'on');
                rawHandle = plot(ax, x, raw(:, component), '.', ...
                    'Color', [0.62 0.64 0.67], 'MarkerSize', 4, ...
                    'DisplayName', 'Observed');
                correctedHandle = plot(ax, x, corrected(:, component), '-', ...
                    'Color', [0.10 0.32 0.60], 'LineWidth', 0.9, ...
                    'DisplayName', 'Offset-corrected');
                obj.addEventLines(ax, detail, component == 1);
                hold(ax, 'off');
                title(ax, sprintf('Detected steps - %s', shortNames{component}));
                ylabel(ax, sprintf('%s (mm)', names{component}));
                if component == 1
                    legend(ax, [rawHandle, correctedHandle], ...
                        'Location', 'best');
                end
                if component == 3
                    xlabel(ax, 'Decimal year');
                end
                obj.applyTimeLimits(ax, x);
            end
        end

        function plotResidual(obj, detail)
            axesList = [obj.App.ResidualsNorthUIAxes, ...
                obj.App.ResidualsEastUIAxes, obj.App.ResidualsUpUIAxes];
            names = {'North', 'East', 'Up'};
            shortNames = {'N', 'E', 'U'};
            x = detail.ts.decimalYear(:);
            if ~isfield(detail.detection, 'finalFit') || ...
                    ~isfield(detail.detection.finalFit, 'residuals')
                for i = 1:3
                    obj.showEmptyAxes(axesList(i), ...
                        'Post-correction residuals are unavailable');
                end
                return;
            end
            residual = 1000 .* detail.detection.finalFit.residuals;
            sigma = obj.residualSigma(detail.detection.finalFit);
            for component = 1:3
                ax = axesList(component);
                obj.resetAxes(ax);
                hold(ax, 'on');
                if isfinite(sigma(component)) && numel(x) > 1
                    patch(ax, [x(1), x(end), x(end), x(1)], ...
                        [-sigma(component), -sigma(component), ...
                        sigma(component), sigma(component)], ...
                        [0.10 0.32 0.60], 'FaceAlpha', 0.08, ...
                        'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
                yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], ...
                    'HandleVisibility', 'off');
                plot(ax, x, residual(:, component), '.', ...
                    'Color', [0.10 0.32 0.60], 'MarkerSize', 4);
                obj.addEventLines(ax, detail, false);
                hold(ax, 'off');
                title(ax, sprintf('Post-correction residual - %s', ...
                    shortNames{component}));
                ylabel(ax, sprintf('%s residual (mm)', names{component}));
                if component == 3
                    xlabel(ax, 'Decimal year');
                end
                obj.applyTimeLimits(ax, x);
            end
        end

        function sigma = residualSigma(~, finalFit)
            sigma = nan(1, 3);
            if isfield(finalFit, 'Sigma') && ...
                    isequal(size(finalFit.Sigma), [3 3]) && ...
                    all(isfinite(finalFit.Sigma), 'all')
                sigma = 1000 .* sqrt(max(real(diag(finalFit.Sigma)), 0)).';
            end
        end

        function plotGlrPower(obj, detail)
            ax = obj.App.GLRPowerUIAxes;
            obj.resetAxes(ax);
            diagnostics = detail.detection.diagnostics;
            x = detail.ts.decimalYear(:);
            if isempty(diagnostics)
                obj.showEmptyAxes(ax, 'No GLR candidate epochs were evaluated');
                return;
            end
            hold(ax, 'on');
            colors = lines(max(1, numel(diagnostics)));
            plotted = 0;
            for i = 1:numel(diagnostics)
                diagnostic = diagnostics(i);
                if ~isfield(diagnostic, 'power') || ...
                        numel(diagnostic.power) ~= numel(x) || ...
                        ~any(isfinite(diagnostic.power))
                    continue;
                end
                plotted = plotted + 1;
                plot(ax, x, diagnostic.power(:), '-', ...
                    'Color', colors(i, :), 'LineWidth', 0.9, ...
                    'DisplayName', sprintf('Iteration %d', i));
                if isfield(diagnostic, 'accepted') && diagnostic.accepted && ...
                        isfinite(diagnostic.candidateIndex)
                    index = round(diagnostic.candidateIndex);
                    if index >= 1 && index <= numel(x)
                        plot(ax, x(index), diagnostic.maxPower, 'o', ...
                            'Color', colors(i, :), ...
                            'MarkerFaceColor', colors(i, :), ...
                            'HandleVisibility', 'off');
                    end
                end
            end
            if isfield(detail.detection, 'threshold') && ...
                    isfinite(detail.detection.threshold)
                yline(ax, detail.detection.threshold, '--', ...
                    sprintf('Threshold %.3f', detail.detection.threshold), ...
                    'Color', [0.78 0.18 0.16], 'LineWidth', 1.0, ...
                    'HandleVisibility', 'off');
            end
            obj.addEventLines(ax, detail, false);
            hold(ax, 'off');
            if plotted == 0
                obj.showEmptyAxes(ax, 'No GLR power curves are available');
                return;
            end
            title(ax, sprintf('GLR Power | %d iteration(s)', plotted));
            xlabel(ax, 'Decimal year');
            ylabel(ax, 'GLR statistic');
            if plotted <= 8
                legend(ax, 'Location', 'best');
            end
            obj.applyTimeLimits(ax, x);
        end

        function plotStepVector(obj, detail)
            ax = obj.App.StepVectorUIAxes;
            obj.resetAxes(ax);
            coefficients = 1000 .* detail.detection.offsetCoefficients;
            if isempty(coefficients)
                obj.showEmptyAxes(ax, 'No accepted step vectors');
                return;
            end
            bars = bar(ax, 1:size(coefficients, 1), coefficients, 'grouped');
            colors = [0.00 0.58 0.74; 0.86 0.37 0.12; 0.18 0.62 0.42];
            for component = 1:min(3, numel(bars))
                bars(component).FaceColor = colors(component, :);
            end
            tableData = detail.detection.offsetTable;
            stdFields = {'NorthStd_mm', 'EastStd_mm', 'UpStd_mm'};
            if istable(tableData) && height(tableData) == size(coefficients, 1) && ...
                    all(ismember(stdFields, tableData.Properties.VariableNames))
                groupCount = size(coefficients, 2);
                groupWidth = min(0.8, groupCount / (groupCount + 1.5));
                hold(ax, 'on');
                for component = 1:groupCount
                    errorX = (1:size(coefficients, 1)) - groupWidth / 2 + ...
                        (2 * component - 1) * groupWidth / (2 * groupCount);
                    errorbar(ax, errorX, coefficients(:, component), ...
                        tableData.(stdFields{component}), ...
                        'LineStyle', 'none', 'Color', [0.28 0.28 0.30], ...
                        'CapSize', 4, 'HandleVisibility', 'off');
                end
                hold(ax, 'off');
            end
            title(ax, 'Step-Vector Components');
            ylabel(ax, 'Signed step (mm)');
            xlabel(ax, 'Accepted offset epoch');
            legend(ax, {'North', 'East', 'Up'}, 'Location', 'best');
            if istable(tableData) && height(tableData) == size(coefficients, 1)
                labels = strings(height(tableData), 1);
                for i = 1:height(tableData)
                    labels(i) = string(tableData.Date(i)) + newline + ...
                        obj.shortEventType(string(tableData.EventClass(i)));
                end
                ax.XTick = 1:numel(labels);
                ax.XTickLabel = cellstr(labels);
                ax.XTickLabelRotation = 25;
            end
        end

        function label = shortEventType(~, eventClass)
            switch eventClass
                case "known-constrained"
                    label = "known";
                case "persistent-candidate"
                    label = "persistent";
                case "low-support-candidate"
                    label = "low support";
                case "transient-pair"
                    label = "transient";
                otherwise
                    label = eventClass;
            end
        end

        function addEventLines(obj, ax, detail, showLabels)
            indices = detail.detection.acceptedIndices(:);
            classes = strings(numel(indices), 1);
            if isfield(detail.detection, 'eventClass') && ...
                    numel(detail.detection.eventClass) == numel(indices)
                classes = string(detail.detection.eventClass(:));
            end
            x = detail.ts.decimalYear(:);
            for i = 1:numel(indices)
                index = indices(i);
                if index < 1 || index > numel(x)
                    continue;
                end
                label = '';
                if showLabels
                    label = sprintf('%d', i);
                end
                line = xline(ax, x(index), ':', label, ...
                    'Color', obj.eventColor(classes(i)), ...
                    'LineWidth', 0.9, 'HandleVisibility', 'off');
                if showLabels
                    line.LabelVerticalAlignment = 'top';
                    line.LabelHorizontalAlignment = 'left';
                end
            end
        end

        function color = eventColor(~, eventClass)
            switch string(eventClass)
                case "known-constrained"
                    color = [0.18 0.45 0.72];
                case "persistent-candidate"
                    color = [0.78 0.18 0.16];
                case "low-support-candidate"
                    color = [0.92 0.55 0.12];
                case "transient-pair"
                    color = [0.55 0.30 0.72];
                otherwise
                    color = [0.45 0.45 0.45];
            end
        end

        function clearPlots(obj)
            axesList = obj.allAxes();
            for i = 1:numel(axesList)
                obj.resetAxes(axesList(i));
            end
            title(obj.App.NorthUIAxes, 'Detected steps - N');
            title(obj.App.EastUIAxes, 'Detected steps - E');
            title(obj.App.UpUIAxes, 'Detected steps - U');
            title(obj.App.ResidualsNorthUIAxes, 'Post-correction residual - N');
            title(obj.App.ResidualsEastUIAxes, 'Post-correction residual - E');
            title(obj.App.ResidualsUpUIAxes, 'Post-correction residual - U');
            title(obj.App.GLRPowerUIAxes, 'GLR Power');
            title(obj.App.StepVectorUIAxes, 'Step-Vector Components');
        end

        function axesList = allAxes(obj)
            axesList = [obj.App.NorthUIAxes, obj.App.EastUIAxes, ...
                obj.App.UpUIAxes, obj.App.ResidualsNorthUIAxes, ...
                obj.App.ResidualsEastUIAxes, obj.App.ResidualsUpUIAxes, ...
                obj.App.GLRPowerUIAxes, obj.App.StepVectorUIAxes];
        end

        function resetAxes(~, ax)
            legend(ax, 'off');
            cla(ax, 'reset');
            grid(ax, 'on');
            box(ax, 'on');
            ax.XTickMode = 'auto';
            ax.YTickMode = 'auto';
            ax.NextPlot = 'replacechildren';
        end

        function showEmptyAxes(obj, ax, message)
            obj.resetAxes(ax);
            text(ax, 0.5, 0.5, message, 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'Color', [0.45 0.45 0.45]);
            ax.XTick = [];
            ax.YTick = [];
        end

        function applyTimeLimits(~, ax, x)
            if numel(x) > 1 && isfinite(x(1)) && isfinite(x(end)) && ...
                    x(end) > x(1)
                xlim(ax, [x(1), x(end)]);
            end
        end

        function updateStationNode(obj, inputFile, status)
            if isempty(obj.Stations)
                return;
            end
            index = find(strcmpi({obj.Stations.inputFile}, inputFile), ...
                1, 'first');
            if isempty(index)
                return;
            end
            node = obj.Stations(index).node;
            if ~isempty(node) && isvalid(node)
                node.Text = sprintf('%s [%s]', ...
                    obj.Stations(index).station, status);
            end
        end

        function station = stationName(~, inputFile)
            [~, stem] = fileparts(inputFile);
            token = regexp(stem, '^(.*)_[012]$', 'tokens', 'once');
            if ~isempty(token)
                stem = token{1};
            end
            station = upper(char(stem));
        end

        function appendLog(obj, message)
            if isempty(message)
                return;
            end
            if ischar(message) || (isstring(message) && isscalar(message))
                messages = {char(string(message))};
            elseif isstring(message)
                messages = cellstr(message(:));
            elseif iscell(message)
                messages = message(:);
            else
                messages = {char(string(message))};
            end
            stamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            for i = 1:numel(messages)
                value = messages{i};
                if ischar(value) || (isstring(value) && isscalar(value))
                    obj.LogLines{end + 1, 1} = sprintf('[%s] %s', ...
                        stamp, char(string(value))); %#ok<AGROW>
                end
            end
            maximumLines = 1000;
            if numel(obj.LogLines) > maximumLines
                obj.LogLines = obj.LogLines(end - maximumLines + 1:end);
            end
            obj.App.TextArea.Value = obj.LogLines;
            drawnow limitrate;
        end

        function setControlsBusy(obj, busy)
            if busy
                state = 'off';
            else
                state = 'on';
            end
            controls = {obj.App.RunallButton, obj.App.RunSelectedButton, ...
                obj.App.ExportFigureButton, ...
                obj.App.DataImportBrowseButton, ...
                obj.App.DataOutputBrowseButton};
            for i = 1:numel(controls)
                control = controls{i};
                if ~isempty(control) && isvalid(control)
                    control.Enable = state;
                end
            end
        end

        function restoreUi(obj, previousTitle)
            try
                if ~isempty(obj.ProgressDialog) && isvalid(obj.ProgressDialog)
                    close(obj.ProgressDialog);
                end
            catch
                % The dialog may already be closed with the app.
            end
            obj.ProgressDialog = [];
            if isempty(obj.App) || ~isvalid(obj.App) || ...
                    isempty(obj.App.GLRUIFigure) || ...
                    ~isvalid(obj.App.GLRUIFigure)
                return;
            end
            obj.IsBusy = false;
            obj.setControlsBusy(false);
            obj.resetStateButton(obj.App.RunallButton);
            obj.resetStateButton(obj.App.RunSelectedButton);
            obj.resetStateButton(obj.App.ExportFigureButton);
            if strcmp(obj.App.GLRUIFigure.Name, 'MATLAB App')
                obj.App.GLRUIFigure.Name = previousTitle;
            end
            drawnow limitrate;
        end

        function resetStateButton(~, button)
            if ~isempty(button) && isvalid(button) && isprop(button, 'Value')
                button.Value = false;
            end
        end
    end
end
