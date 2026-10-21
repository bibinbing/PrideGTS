function results = StationsDisError(data,ref)
% function: Calculate the distance correlation between the reference
%           station and each measuring station
% Parameters
% data: The filtered data matrix,is .mat, row is time coordinate
%       the column is station.
% ref : Number of reference station column.
% Usage: results = StationsDisError(data,1).

data_size = size(data);
rows = data_size(2);

progressBar = waitbar(0,'Processing...');

for i = 1:rows
    results(i) = distcorr(data(:,ref),data(:,i));
    waitbar(i / rows, progressBar);
end

close(progressBar);
end