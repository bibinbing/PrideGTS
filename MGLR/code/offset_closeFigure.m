function offset_closeFigure(fig)
% Closes a figure only when its graphics handle is still valid.

if ~isempty(fig) && isgraphics(fig, 'figure')
    close(fig);
end

end
