function VBICAResult = RunvbICA(X,NumReconComp,NumGuassModel,Source_type,SEE_PDFS,tol,max_steps,ARD)

VBICAResult = struct;

[m,n] = size(X);
states = repmat(NumGuassModel,1,n);
X1 = X'; % 13*1553

net = vbica2(X1(:,1:m),states,Source_type,SEE_PDFS,tol,max_steps,1,ARD);

numOfIC = n;
Contribution = zeros(numOfIC, 1);
for i = 1:numOfIC
    Contribution(i) = (var(net.recon(i,:),0,2)/sum(var(net.recon, 0, 2))) * 100; % Total variance across features
end

[SortContri, Idx] = sort(Contribution, 'descend');
A = net.mixmatrix(:, Idx);
S = net.recon(Idx, :);

cme = A(:,1:NumReconComp) * S(1:NumReconComp,:);
FilteredData = X - cme';

VBICAResult.Contribution = SortContri;
VBICAResult.FilteredData = FilteredData;
VBICAResult.OriginTimeSeries = X;
VBICAResult.SpatialResponsesMatrix = A; % ATTENTION!!!
VBICAResult.Components = S';
end