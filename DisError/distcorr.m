function dcor = distcorr(x, y)
%DISTCOR Biased sample distance correlation using base MATLAB only.
% Rows are observations; columns are variables. Remove rows with NaN in
% either input, as in the original PrideGTS implementation.
% The estimator and double centering follow the original distcorr.m by
% Shen Liu (18 Jan 2013); Euclidean distances no longer require pdist2.
% Constant inputs or fewer than two paired observations return NaN because
% their sample distance correlation has a zero/undefined denominator.

validateattributes(x, {'numeric'}, {'real', '2d'}, mfilename, 'x');
validateattributes(y, {'numeric'}, {'real', '2d'}, mfilename, 'y');
if size(x, 1) ~= size(y, 1)
    error('PrideGTS:distcorr:RowMismatch', ...
        'Inputs must have the same number of rows.');
end
missing = any(isnan(x), 2) | any(isnan(y), 2);
x = double(x(~missing, :));
y = double(y(~missing, :));
if any(isinf(x(:))) || any(isinf(y(:)))
    error('PrideGTS:distcorr:InfiniteInput', 'Inputs must not contain Inf.');
end
if size(x, 1) < 2 || isempty(x) || isempty(y)
    dcor = NaN;
    return
end
A = centeredDistances(x);
B = centeredDistances(y);
% The common factor 1/N^2 cancels in the correlation ratio.
varx = sum(A(:).^2);
vary = sum(B(:).^2);
if varx == 0 || vary == 0
    dcor = NaN;
    return
end
covxy = sum(A(:).*B(:));
ratio = (covxy / sqrt(varx)) / sqrt(vary);
dcor = sqrt(min(1, max(0, ratio)));
end

function A = centeredDistances(x)
n = size(x, 1);
D = zeros(n, n);
for k = 1:size(x, 2)
    % Direct differences avoid cancellation for large coordinate origins.
    D = hypot(D, bsxfun(@minus, x(:, k), x(:, k).'));
end
rowMean = mean(D, 2);
colMean = mean(D, 1);
A = bsxfun(@minus, bsxfun(@minus, D, rowMean), colMean) ...
    + mean(colMean);
end
