function outFile = offset_plotPower(ts, detection, outFile, opts)
% Plots the score used to locate possible offsets.

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 520]);
figureCleanup = onCleanup(@() offset_closeFigure(fig));
x = ts.decimalYear;

if isempty(detection.diagnostics)
    text(0.5, 0.5, 'No candidate epochs evaluated', 'HorizontalAlignment', 'center');
    axis off;
else
    hold on;
    diagnostics = detection.diagnostics;
    progression = linspace(0, 1, max(numel(diagnostics), 1)).';
    colors = [0.72 - 0.40 .* progression, ...
              0.60 - 0.38 .* progression, ...
              0.82 - 0.16 .* progression];
    for k = 1:numel(diagnostics)
        if ~isfield(diagnostics(k), 'power') || isempty(diagnostics(k).power)
            continue;
        end
        plot(x, diagnostics(k).power, '-', ...
            'Color', colors(k, :), 'LineWidth', 0.9);
        if diagnostics(k).accepted && isfinite(diagnostics(k).candidateIndex)
            idx = diagnostics(k).candidateIndex;
            scatter(x(idx), diagnostics(k).maxPower, 30, ...
                colors(k, :), 'filled', 'MarkerEdgeColor', 'k');
        end
    end
    yline(detection.pointwiseThreshold, '--', ...
        'Color', [0.55 0.57 0.60], 'LineWidth', 0.9);
    yline(detection.threshold, '-', ...
        'Color', [0.93 0.53 0.12], 'LineWidth', 1.2);
    offset_addLines(x, detection.acceptedIndices, ...
        detection.applyCorrection, detection.eventClass);
    xlabel('Decimal year');
    ylabel('Multivariate GLR statistic');
    title(sprintf(['%s station-level GLR | %d iterations | pointwise ', ...
        '%.3g | station %.3g (alpha %.4g)'], ...
        ts.station, numel(diagnostics), detection.pointwiseThreshold, ...
        detection.threshold, detection.alpha), ...
        'Interpreter', 'none');
    grid on;
end

offset_savepng(fig, outFile, opts);
clear figureCleanup

end
