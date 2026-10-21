function h = offset_addLines(x, acceptedIndices, applyCorrection, eventClass)
% Adds classified vertical markers at detected candidate epochs.

h = gobjects(0);
if isempty(acceptedIndices)
    return;
end
if nargin < 3 || isempty(applyCorrection)
    applyCorrection = true(numel(acceptedIndices), 1);
end
if nargin < 4 || isempty(eventClass)
    eventClass = repmat("persistent-candidate", numel(acceptedIndices), 1);
end
applyCorrection = logical(applyCorrection(:));
eventClass = string(eventClass(:));
if numel(applyCorrection) ~= numel(acceptedIndices) || ...
        numel(eventClass) ~= numel(acceptedIndices)
    error('offset:EventLineAlignment', ...
        'Candidate epochs, correction flags, and event classes must align.');
end

yl = ylim;
h = gobjects(numel(acceptedIndices), 1);
for i = 1:numel(acceptedIndices)
    xv = x(acceptedIndices(i));
    if eventClass(i) == "known-constrained"
        color = [0.43 0.24 0.70];
        style = '-.';
        width = 1.2;
    elseif applyCorrection(i)
        color = [0.78 0.18 0.16];
        style = '-';
        width = 1.0;
    elseif eventClass(i) == "transient-pair"
        color = [0.93 0.53 0.12];
        style = '--';
        width = 1.0;
    else
        color = [0.62 0.64 0.67];
        style = ':';
        width = 0.8;
    end
    h(i) = line([xv xv], yl, 'Color', color, ...
        'LineStyle', style, 'LineWidth', width);
end
ylim(yl);

end
