function value = offset_chi2crit(alpha, dof)
% Computes the upper-tail chi-square threshold for the GLR statistic.

if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0 || alpha >= 1
    error('offset:InvalidAlpha', ...
        'The GLR significance level alpha must lie strictly between 0 and 1.');
end
if ~isscalar(dof) || ~isfinite(dof) || dof <= 0 || dof ~= fix(dof)
    error('offset:InvalidChiSquareDof', ...
        'The chi-square degrees of freedom must be a positive integer.');
end

if exist('chi2inv', 'file') == 2
    value = chi2inv(1 - alpha, dof);
else
    value = 2 .* gammaincinv(1 - alpha, dof ./ 2);
end

end
