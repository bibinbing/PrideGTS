function name = offset_safe(value)
% Converts text to a safe file-name stem.

name = regexprep(char(value), '[^A-Za-z0-9_.-]', '_');
if isempty(name)
    name = 'station';
end

end
