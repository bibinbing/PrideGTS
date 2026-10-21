function mjd = yymmmdd2mjd(dateStr)
%YYMMMDD2MJD Convert an NGL steps.txt date token (e.g. '17MAY04') to MJD.
%
%   mjd = YYMMMDD2MJD(dateStr)
%
%   NGL's steps.txt (https://geodesy.unr.edu/NGLStationPages/steps.txt)
%   encodes each step's date as a 2-digit-year + 3-letter-month + 2-digit
%   day token, e.g. '17MAY04' = 2017-05-04. Two-digit years are resolved
%   assuming the GNSS era (1970-2069): yy >= 70 -> 19yy, otherwise 20yy.
%   This is unambiguous for all dates found in real GNSS steps catalogues
%   (1994-present).
%
%   Example:
%     yymmmdd2mjd('17MAY04')   % -> 57877

dateStr = strtrim(dateStr);
if numel(dateStr) ~= 7
    error('yymmmdd2mjd:badFormat', 'Expected a 7-character token like ''17MAY04'', got ''%s''.', dateStr);
end

yy  = str2double(dateStr(1:2));
mon = upper(dateStr(3:5));
dd  = str2double(dateStr(6:7));

if yy >= 70
    yyyy = 1900 + yy;
else
    yyyy = 2000 + yy;
end

months = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
mm = find(strcmp(months, mon), 1);
if isempty(mm)
    error('yymmmdd2mjd:badMonth', 'Unrecognized month token in ''%s''.', dateStr);
end

mjd = ymd2mjd(yyyy, mm, dd);

end
