function lines = offset_txt(fileName)
% Reads a text file into a list of lines.

text = fileread(fileName);
lines = regexp(text, '\r\n|\n|\r', 'split')';
if ~isempty(lines) && isempty(lines{end})
    lines(end) = [];
end

end
