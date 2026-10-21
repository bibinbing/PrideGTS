% This demonstrates vbICA2 on 3 dimensional data produced by mixing 2
% sources and adding 5% Gaussian noise.

load test.mat
data = alldata(:,[501:1000,3001:4500]);

fprintf('Press any key to plot 2000-vector data distribution...')
pause
fprintf('\n')
figure(1)
plot3(data(1,:),data(2,:),data(3,:),'.')

fprintf('ICA network will learn 3 MoG sources with 3 states each, \n')
fprintf('and with isotropic noise on 1000 data points.\n')
fprintf('Press any key to continue...')
pause
fprintf('\n')
figure(2)
net = vbica2(data(:,1:1000),[3 3 3],'g',1);


fprintf('Learning has finished. Press any key to plot reconstructions...')
pause
fprintf('\n')
figure(3)
plotsignals(net.recon)
fprintf('One of the sources has been suppressed by ARD - two \n')
fprintf('significant sources found.\n')

fprintf('Press any key to unmix rest of the data...')
pause
fprintf('\n')
recon = unmix2(data,net);


fprintf('Press any key to plot underlying distribution...')
pause
fprintf('\n')
figure(4)
plot(recon.signals(1,:),recon.signals(2,:),'.')

fprintf('Press any key to plot (significant) source models...')
pause
fprintf('\n')
figure(5)
subplot(2,2,1),hist(sources{2}(1,:),100)
subplot(2,2,2),hist(sources{2}(2,:),100)
subplot(2,2,3),plotMoG(net.sources{1})
subplot(2,2,4),plotMoG(net.sources{2})
fprintf('Sources may be permuted and scaled...\n\n')
fprintf('The End\n')