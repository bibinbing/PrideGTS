function outFile = offset_plotCompare(ts, detection, correctedNeu, outFile, opts)
% Plots raw and corrected NEU values for comparison.

x = ts.decimalYear;
raw = ts.neu .* 1000;
corr = correctedNeu .* 1000;
labels = {'North (mm)', 'East (mm)', 'Up (mm)'};
componentColors = [0.00 0.58 0.74; ...
                   0.86 0.37 0.12; ...
                   0.18 0.62 0.42];
beforeColor = [0.62 0.64 0.67];
afterColor = [0.10 0.32 0.60];

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 820]);
figureCleanup = onCleanup(@() offset_closeFigure(fig));
for k = 1:3
    subplot(3, 1, k);
    plot(x, raw(:, k), '.', 'Color', beforeColor, 'MarkerSize', 5); hold on;
    plot(x, corr(:, k), '-', 'Color', afterColor, 'LineWidth', 0.9);
    offset_addLines(x, detection.acceptedIndices, ...
        detection.applyCorrection, detection.eventClass);
    yLabel = ylabel(labels{k});
    yLabel.Color = componentColors(k, :);
    grid on;
    if k == 1
        title(sprintf('%s: raw (gray) vs corrected (blue)', ...
            ts.station), 'Interpreter', 'none');
    end
    if k == 3
        xlabel('Decimal year');
    end
end
offset_savepng(fig, outFile, opts);
clear figureCleanup

end
