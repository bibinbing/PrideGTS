function mjd = ymd2mjd(yyyy, mm, dd)
%YMD2MJD Convert a Gregorian calendar date to a Modified Julian Date.
%
%   mjd = YMD2MJD(yyyy, mm, dd)
%
%   Dependency-free implementation of the Fliegel & Van Flandern
%   algorithm. Works in base MATLAB with no toolbox requirement, so it
%   is safe to use regardless of which MATLAB edition PrideGTS is run
%   with.
%
%   Example:
%     ymd2mjd(2000,1,1)   % -> 51544

a = floor((14 - mm) / 12);
y = yyyy + 4800 - a;
m = mm + 12*a - 3;

jdn = dd + floor((153*m + 2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045;
mjd = jdn - 2400001;

end
