function mjd = decyear2mjd(decYear)
%DECYEAR2MJD Convert a CEA-style decimal year to a Modified Julian Date.
%
%   mjd = DECYEAR2MJD(decYear)
%
%   CEA (First Crust Monitoring and Application Center) time series
%   files encode epochs as "YYYY.DECM" using the convention
%
%       decYear = year + (DOY - 0.5) / daysInYear(year)
%
%   i.e. each day's fractional-year value is referenced to local noon.
%   This matches the "YYYY.DECM" column found in CEA coordinate files
%   exactly (verified against BJFS.txt: 2011.19041 <-> DOY 70 in 2011;
%   2015.90274 <-> DOY 330 in 2015).
%
%   Example:
%     decyear2mjd(2011.19041)   % -> 55631

yyyy = floor(decYear);
frac = decYear - yyyy;

if (mod(yyyy,4)==0 && mod(yyyy,100)~=0) || mod(yyyy,400)==0
    daysInYear = 366;
else
    daysInYear = 365;
end

doy = round(frac * daysInYear + 0.5);

mjd0 = ymd2mjd(yyyy, 1, 1);   % MJD of Jan 1st of that year
mjd  = mjd0 + (doy - 1);

end
