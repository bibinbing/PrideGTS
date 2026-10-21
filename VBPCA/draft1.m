clc
clear
%%
load('C:\Users\26403\Desktop\GmisFormatData.mat')
data = GmisStruct.dn;
OUT = vbpca_full(data, 10, 'display', 1, 'bias', 0, 'niter_broadprior', 100);

%% comtri
row_var = OUT.Va;
total_var = sum(row_var);
%This is contribution ratios.
variance_ratios = (row_var / total_var).*100;

% Sort value and index.
[variance_ratios,ContriIdex] = sort(variance_ratios,"descend");
total_ratios =cumsum(variance_ratios);

figure;
hold on;
plot(variance_ratios, 'k-', 'LineWidth', 1);
plot(total_ratios, 'r-', 'LineWidth', 1);
xlabel('Principal Component Index');
ylabel('Variance Ratio (%)');
title('Comtribution Ratios of Principal Components');
grid on;
plot(1:length(total_ratios), total_ratios, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
legend({'Individual Variance Ratio', 'Cumulative Variance Ratio'}, 'Location', 'best');
hold off;

%%
A = OUT.A(:,ContriIdex);
S = OUT.S(ContriIdex,:); 
Xc = data ;
CME   = A(:,1)*S(1,:);         
filtX = Xc - CME - OUT.Mu ;           

down_ratio = zeros(size(Xc,2), 1);
rmse_before = zeros(size(Xc,2), 1);
rmse_after = zeros(size(Xc,2), 1);
for i = 1:size(Xc,2)
    index=find(~isnan(Xc(:,i)));
    rmse_before(i) = rmse(Xc(index,i));          % Before RMSE
    rmse_after(i) = rmse(filtX(index,i));    % After RMSE
    down_ratio(i,1) = ((rmse_before(i) - rmse_after(i)) / rmse_before(i)) * 100;
end

figure
hold on
ylabel("mm");
bar(rmse_before,'red');
bar(rmse_after,'blue');
plot(down_ratio./100,'r-','LineWidth',2,'DisplayName','ratios');
legend("滤波前","滤波后","下降率");
grid on
hold off
