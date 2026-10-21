function [pxx, frequency] = pridegts_lomb_psd(values, sampleTime)
%PRIDEGTS_LOMB_PSD Single-channel Lomb-Scargle PSD in base MATLAB.
% This implements only the two-vector call used by plotResidualSpectra_.
% It is not a replacement for every syntax of Signal Processing plomb.
% Time is numeric in any consistent unit; outputs use cycles/time-unit and
% values^2/(cycles/time-unit). For years and mm, PSD units are mm^2 yr.
% NaN pairs are omitted before grid AND normalization are calculated, as in
% PrideGTS's caller. This need not reproduce plomb(x,t) with NaNs still in x.
% Remaining times must be strictly increasing.
% The caller retains responsibility for linear detrending.
%
% Grid: dt=(t(end)-t(1))/(N-1), df=1/(4*N*dt), f=(1:2*N)*df.
% Classical mean-subtracted Lomb-Scargle projection, one-sided PSD scaling.
% Reference: Scargle (1982), ApJ 263:835-853.
% Public interface/formula: https://www.mathworks.com/help/signal/ref/plomb.html

validateattributes(values, {'numeric'}, {'real', 'vector'}, mfilename, 'values');
validateattributes(sampleTime, {'numeric'}, {'real', 'vector'}, mfilename, 'sampleTime');
if numel(values) ~= numel(sampleTime)
    error('PrideGTS:lomb:SizeMismatch', 'Values and times must have equal lengths.');
end
x = double(values(:));
t = double(sampleTime(:));
keep = ~isnan(x) & ~isnan(t);
x = x(keep);
t = t(keep);
if numel(t) < 3
    error('PrideGTS:lomb:TooFewSamples', 'At least three valid samples are required.');
end
if any(~isfinite(x)) || any(~isfinite(t)) || any(diff(t) <= 0)
    error('PrideGTS:lomb:InvalidSamples', ...
        'Samples must be finite and times strictly increasing.');
end
n = numel(t);
t = t - t(1);
dt = t(end)/(n - 1);
frequency = (1:2*n).'/(4*n*dt);
x = x - mean(x);
pxx = zeros(size(frequency));
% Bound the temporary trigonometric matrices for decade-long daily series.
blockSize = 128;
for first = 1:blockSize:numel(frequency)
    idx = first:min(first + blockSize - 1, numel(frequency));
    phase = 2*pi*(t*frequency(idx).');
    shift = 0.5*atan2(sum(sin(2*phase), 1), sum(cos(2*phase), 1));
    phase = bsxfun(@minus, phase, shift);
    C = cos(phase);
    S = sin(phase);
    cc = sum(C.^2, 1);
    ss = sum(S.^2, 1);
    xc = x.'*C;
    xs = x.'*S;
    powerC = zeros(size(cc));
    powerS = zeros(size(ss));
    % At exact regular-sampling Nyquist, one basis can have zero rank.
    goodC = cc > n*eps;
    goodS = ss > n*eps;
    powerC(goodC) = xc(goodC).^2 ./ cc(goodC);
    powerS(goodS) = xs(goodS).^2 ./ ss(goodS);
    pxx(idx) = dt*(powerC + powerS).';
end
end
