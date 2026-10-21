function offset_savepng(fig, outFile, opts)
% Exports one PNG atomically and releases its graphics resources.

figureCleanup = onCleanup(@() offset_closeFigure(fig));
[folder, ~, extension] = fileparts(outFile);
if isempty(folder)
    folder = pwd;
end
if isempty(extension)
    extension = '.png';
end
temporaryFile = [tempname(folder), extension];
fileCleanup = onCleanup(@() removeTemporaryFile(temporaryFile));

drawnow;
exportgraphics(fig, temporaryFile, ...
    'Resolution', round(opts.figures.resolution), ...
    'BackgroundColor', 'white');
[moved, message] = movefile(temporaryFile, outFile, 'f');
if ~moved
    error('offset:FigureMoveFailed', ...
        'Could not finalize figure %s: %s', outFile, message);
end

offset_closeFigure(fig);
drawnow;
clear figureCleanup fileCleanup

end

function removeTemporaryFile(fileName)
if isfile(fileName)
    delete(fileName);
end

end
