function outFile = offset_write(ts, correctedNeu, opts, detection)
% Writes the MOM station output.

if nargin < 4
    detection = [];
end

if ~isfield(ts, 'sourceFormat') || ~strcmp(ts.sourceFormat, 'MOM')
    error('offset:UnsupportedWriteFormat', 'Only MOM output is supported.');
end

outFile = offset_writeMOM(ts, correctedNeu, detection, opts);

end
