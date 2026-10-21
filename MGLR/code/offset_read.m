function ts = offset_read(fileName)
% Reads one three-component MOM station.

[~, ~, extension] = fileparts(fileName);
if ~strcmpi(extension, '.mom')
    error('offset:MOMExtension', 'Only .mom input files are supported: %s.', fileName);
end

ts = offset_readMOM(fileName);
[ts, ts.quality] = offset_cleanSeries(ts);
ts.count = numel(ts.mjd);

end
