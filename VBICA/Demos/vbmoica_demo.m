% This demonstrates vbMoICA on 3 dimensional data containing 3 clusters
% of 3 different intrinsic dimensionalities.


load test.mat
addpath ../vbMoICA


fprintf('Press any key to plot 6000-vector data distribution...')
pause
fprintf('\n')
figure(1)
plot3(alldata(1,[1:500,1501:3000]),alldata(2,[1:500,1501:3000]), ...
      alldata(3,[1:500,1501:3000]),'b.')
hold on
plot3(alldata(1,[501:1000,3001:4500]),alldata(2,[501:1000,3001:4500]), ...
      alldata(3,[501:1000,3001:4500]),'r.')
plot3(alldata(1,[1001:1500,4501:6000]),alldata(2,[1001:1500,4501:6000]), ...
      alldata(3,[1001:1500,4501:6000]),'g.')

fprintf('vbMoICA network will learn 3 ICA networks each with 3 MoG \n')
fprintf('with 3 MoG states each, \n')
fprintf('and with isotropic noise on 1500 data points.\n')
fprintf('Press any key to continue...')
pause
fprintf('\n')
figure(2)
mixnet = vbmoica(alldata(:,1:500),states(1:3));


fprintf('Learning has finished. Press any key to plot reconstructions...')
pause
fprintf('\n')
figure(3)
ica = mixnet.nets;
plotsignals([ica{1}.recon;ica{2}.recon;ica{3}.recon],3)
fprintf('Correct intrinsic dimensionalities found.\n')

fprintf('Press any key to unmix rest of the data...')
pause
fprintf('\n')
[recons eta] = unmix(alldata,mixnet);

fprintf('Press any key to plot partition...')
pause
fprintf('\n')
clrs = ['b.','r.','g.'];
figure(1),subplot(2,1,1)
plot3(alldata(1,[1:500,1501:3000]),alldata(2,[1:500,1501:3000]), ...
      alldata(3,[1:500,1501:3000]),'b.')
hold on
plot3(alldata(1,[501:1000,3001:4500]),alldata(2,[501:1000,3001:4500]), ...
      alldata(3,[501:1000,3001:4500]),'r.')
plot3(alldata(1,[1001:1500,4501:6000]),alldata(2,[1001:1500,4501:6000]), ...
      alldata(3,[1001:1500,4501:6000]),'g.')
subplot(2,1,2)
[mx,a] = max(eta);
j = find(a==1);
plot3(alldata(1,j),alldata(2,j),alldata(3,j),'b.')
hold on
j = find(a==2);
plot3(alldata(1,j),alldata(2,j),alldata(3,j),'r.')
hold on
j = find(a==3);
plot3(alldata(1,j),alldata(2,j),alldata(3,j),'g.')


fprintf('Press any key to plot all source models...')
pause
fprintf('\n')
figure(4)
subplot(3,3,1),hist(sources{1}(1,:),100)
subplot(3,3,2),hist(sources{1}(2,:),100)
subplot(3,3,3),hist(sources{1}(3,:),100)
subplot(3,3,4),hist(sources{2}(1,:),100)
subplot(3,3,5),hist(sources{2}(2,:),100)
subplot(3,3,7),hist(sources{3}(1,:),100)
figure(5)
steps = 0;
for i=1:3
  for j=1:3
    steps = steps+1;
    subplot(3,3,steps),plotMoG(ica{i}.sources{j})
  end
end

fprintf('Sources may be permuted and scaled...\n\n')
fprintf('The End\n')