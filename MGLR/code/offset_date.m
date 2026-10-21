function t = offset_date(mjd)
% Converts Modified Julian Date values to datetime values.

t = datetime(1858, 11, 17, 0, 0, 0) + days(mjd);
t.Format = 'yyyy-MM-dd';

end
