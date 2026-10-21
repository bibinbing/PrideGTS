%% Variational Matrix Low-Rank Matrix Completion
clear
clc

%% Time series data matrix
load("GmisFormatData.mat");
Y = GmisStruct.du;
Y(isnan(Y)) = 0;
P = ~isnan(Y);

%% Run VBMC

fprintf('Running VB matrix completion...\n');
tic
[X_hat, A_hat, B_hat,t,gamma] = VBMC(P,Y);
t_total = toc;
