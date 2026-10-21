function outFile = offset_plotVectors(ts, detection, outFile, opts)
% Plots the detected NEU offset components.

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 520]);
figureCleanup = onCleanup(@() offset_closeFigure(fig));
coeffMm = detection.offsetCoefficients .* 1000;

if isempty(coeffMm)
    text(0.5, 0.5, 'No significant offsets detected', 'HorizontalAlignment', 'center');
    axis off;
else
    hBars = bar(1:size(coeffMm, 1), coeffMm, 'grouped');
    componentColors = [0.00 0.58 0.74; ...
                       0.86 0.37 0.12; ...
                       0.18 0.62 0.42];
    for k = 1:min(numel(hBars), 3)
        hBars(k).FaceColor = componentColors(k, :);
    end
    hold on;
    if isfield(detection, 'offsetTable') && ...
            all(ismember({'NorthStd_mm', 'EastStd_mm', 'UpStd_mm'}, ...
            detection.offsetTable.Properties.VariableNames))
        standardErrors = detection.offsetTable{:, ...
            {'NorthStd_mm', 'EastStd_mm', 'UpStd_mm'}};
        for k = 1:min(numel(hBars), 3)
            errorbar(hBars(k).XEndPoints, coeffMm(:, k), standardErrors(:, k), ...
                'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);
        end
    end
    yline(0, '-', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
    grid on;
    ylabel('Signed offset component (mm)');
    xlabel('Detected offset epoch');
    title(sprintf( ...
        '%s events (%d known, %d GLR, %d corrected, %d transient pairs)', ...
        ts.station, detection.knownConstraintCount, detection.candidateCount, ...
        detection.correctedOffsetCount, ...
        detection.transientPairCount), 'Interpreter', 'none');
    legend(hBars, {'North', 'East', 'Up'}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
    dateValues = ts.time(detection.acceptedIndices);
    dateValues.Format = 'yyyy-MM-dd';
    dates = string(dateValues);
    set(gca, 'XTick', 1:numel(dates), 'XTickLabel', dates);
    xtickangle(30);
end

offset_savepng(fig, outFile, opts);
clear figureCleanup

end
