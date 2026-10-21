function outFile = offset_plotResiduals(ts, detection, outFile, opts)
% Plots residuals before and after adding detected offsets.

x = ts.decimalYear;
rawResiduals = detection.baseFit.residuals .* 1000;
finalResiduals = detection.finalFit.residuals .* 1000;
labels = {'North residual (mm)', 'East residual (mm)', 'Up residual (mm)'};
componentColors = [0.00 0.58 0.74; ...
                   0.86 0.37 0.12; ...
                   0.18 0.62 0.42];
beforeColor = [0.62 0.64 0.67];
afterColor = [0.10 0.32 0.60];
sigmaMm = nan(1, 3);
if isfield(detection, 'finalFit') && ...
        isfield(detection.finalFit, 'Sigma') && ...
        isequal(size(detection.finalFit.Sigma), [3 3]) && ...
        all(isfinite(detection.finalFit.Sigma), 'all')
    sigmaMm = 1000 .* sqrt(max(real(diag(detection.finalFit.Sigma)), 0)).';
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 820]);
figureCleanup = onCleanup(@() offset_closeFigure(fig));
for k = 1:3
    subplot(3, 1, k);
    hold on;
    if isfinite(sigmaMm(k)) && sigmaMm(k) > 0 && numel(x) > 1
        patch([x(1), x(end), x(end), x(1)], ...
            [-sigmaMm(k), -sigmaMm(k), sigmaMm(k), sigmaMm(k)], ...
            afterColor, 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    yline(0, '--', 'Color', [0.32 0.33 0.35], ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(x, rawResiduals(:, k), '.', 'Color', beforeColor, 'MarkerSize', 5);
    plot(x, finalResiduals(:, k), '-', 'Color', afterColor, 'LineWidth', 0.9);
    offset_addLines(x, detection.acceptedIndices, ...
        detection.applyCorrection, detection.eventClass);
    yLabel = ylabel(labels{k});
    yLabel.Color = componentColors(k, :);
    grid on;
    if k == 1
        title(sprintf('%s residuals: base (gray), final (blue), shaded +/-1 sigma', ...
            ts.station), 'Interpreter', 'none');
    end
    if k == 3
        xlabel('Decimal year');
    end
end
offset_savepng(fig, outFile, opts);
clear figureCleanup

end
