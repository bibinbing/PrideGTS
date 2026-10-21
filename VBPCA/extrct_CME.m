%% comtri
row_var = sum(OUT.S.^2,2);
total_var = sum(row_var);
variance_ratios = (row_var / total_var).*100;
total_ratios =cumsum(variance_ratios);

figure;
hold on;
plot(variance_ratios, 'k-', 'LineWidth', 1);
plot(total_ratios, 'r-', 'LineWidth', 1);
xlabel('Principal Component Index');
ylabel('Variance Ratio (%)');
title('Low Rank Matrix Variance Ratios of Principal Components');
grid on;
plot(1:length(total_ratios), total_ratios, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); 
legend({'Individual Variance Ratio', 'Cumulative Variance Ratio'}, 'Location', 'best');
hold off;

