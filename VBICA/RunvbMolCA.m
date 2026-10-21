clc
clear

load('C:\Users\26403\Desktop\temp\GmisFormatData.mat')
alldata = GmisStruct.dn;
[m,n] = size(alldata);
data = alldata'; % 13*1553


% 设置参数
states = {[3 3]};


% 运行 vbMoICA
mixnet = vbmoica(data(:,1:m), states);

%%
numOfIC = size(mixnet.nets{1, 1}.mixmatrix ,2);
Contribution = zeros(numOfIC, 1);
XX = mixnet.nets{1, 1}.mixmatrix * mixnet.nets{1, 1}.recon;

for i = 1:numOfIC
    Contribution(i) = (var(XX(i,:),0,2)/sum(var(XX, 0, 2))) * 100; % Total variance across features
end
