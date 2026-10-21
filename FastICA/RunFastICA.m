function FastICAResult = RunFastICA(X, NumReconComp)

FastICAResult = struct;

% Input validation
if isempty(X)
    error('Input data X cannot be empty.');
end


try
    [S, A, W] = FastICA(X');
catch
    warning('FastICA failed.');
    return
end

numOfIC = size(S, 1);
Contribution = zeros(numOfIC, 1);
XX = A' * S; % Project
for i = 1:numOfIC
    Contribution(i) = (var(XX(i,:),0,2)/sum(var(XX, 0, 2))) * 100; % Total variance across features
end

[SortContri, Idx] = sort(Contribution, 'descend');
A = A(:, Idx);
S = S(Idx, :);

cme = A(:,1:NumReconComp) * S(1:NumReconComp,:);
FilteredData = X - cme';
FastICAResult.Contribution = SortContri;
FastICAResult.FilteredData = FilteredData;
FastICAResult.OriginTimeSeries = X;
FastICAResult.SpatialResponsesMatrix = A; % ATTENTION!!!
FastICAResult.Components = S';

end
