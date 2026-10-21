function decimalYear = offset_decyr(t)
% Converts datetime values to decimal years.

y = year(t);
startOfYear = datetime(y, 1, 1);
startOfNextYear = datetime(y + 1, 1, 1);
decimalYear = y + days(t - startOfYear) ./ days(startOfNextYear - startOfYear);
decimalYear = decimalYear(:);

end
