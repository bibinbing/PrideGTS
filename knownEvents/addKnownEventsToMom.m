function addKnownEventsToMom(station, momDir, outDir, varargin)
%ADDKNOWNEVENTSTOMOM Import known offset epochs into a station's mom files.
%
%   addKnownEventsToMom(station, momDir, outDir, 'StepsFile', f1)
%   addKnownEventsToMom(station, momDir, outDir, 'CEAFile', f2)
%   addKnownEventsToMom(station, momDir, outDir, 'StepsFile', f1, 'CEAFile', f2)
%
%   Combines known-offset epochs from an NGL steps.txt (global stations,
%   incl. Hong Kong) and/or a CEA header-annotated coordinate file
%   (Mainland China stations, which are absent from NGL) and prepends
%   them as a header comment block to the three component .mom files:
%
%       <momDir>/<station>_0.mom   (N)
%       <momDir>/<station>_1.mom   (E)
%       <momDir>/<station>_2.mom   (U)
%
%   NGL-sourced events are not component-specific (an equipment change or
%   nearby earthquake potentially offsets all three components), so they
%   are written to all three files. CEA-sourced events are annotated per
%   component in the source file, so they are only written to the
%   matching component's file.
%
%   Outputs are written to outDir; momDir is left untouched.
%
%   Example:
%     addKnownEventsToMom('HKSC', 'mom_raw', 'mom_annotated', ...
%         'StepsFile', 'steps.txt');

p = inputParser;
p.addParameter('StepsFile', '');
p.addParameter('CEAFile', '');
p.parse(varargin{:});
stepsFile = p.Results.StepsFile;
ceaFile   = p.Results.CEAFile;

if isempty(stepsFile) && isempty(ceaFile)
    error('addKnownEventsToMom:noSource', ...
        'Provide at least one of ''StepsFile'' (NGL) or ''CEAFile'' (CEA).');
end

if ~isfolder(outDir)
    mkdir(outDir);
end

compSuffix = {'0','1','2'};
compName   = {'N','E','U'};

% --- gather events -----------------------------------------------------
ngEvents = table();
if ~isempty(stepsFile)
    ngEvents = parseNGLSteps(stepsFile, station);
    fprintf('[NGL]  %d known step(s) found for %s in %s\n', ...
        height(ngEvents), station, stepsFile);
end

ceaEvents = table();
if ~isempty(ceaFile)
    ceaEvents = parseCEAHeader(ceaFile);
    fprintf('[CEA]  %d known offset(s) found in header of %s\n', ...
        height(ceaEvents), ceaFile);
end

if height(ngEvents) == 0 && height(ceaEvents) == 0
    warning('addKnownEventsToMom:noEvents', ...
        'No known events found for station %s.', station);
end

% --- write each component -----------------------------------------------
for c = 1:3
    momFile = fullfile(momDir, sprintf('%s_%s.mom', station, compSuffix{c}));
    if ~isfile(momFile)
        warning('addKnownEventsToMom:missingComponent', ...
            'Missing %s, skipped.', momFile);
        continue;
    end

    thisEvents = table('Size', [0 3], ...
        'VariableTypes', {'double','string','string'}, ...
        'VariableNames', {'MJD','Source','Description'});

    if height(ngEvents) > 0
        t = table(ngEvents.MJD, repmat("NGL", height(ngEvents), 1), ...
            strcat(ngEvents.Type, ": ", ngEvents.Description), ...
            'VariableNames', {'MJD','Source','Description'});
        thisEvents = [thisEvents; t]; %#ok<AGROW>
    end

    if height(ceaEvents) > 0
        sel = ceaEvents.Component == compName{c};
        if any(sel)
            sub = ceaEvents(sel, :);
            t = table(sub.MJD, repmat("CEA", height(sub), 1), ...
                strcat("offset ", string(sub.Offset_m), " +/- ", string(sub.Sigma_m), " m"), ...
                'VariableNames', {'MJD','Source','Description'});
            thisEvents = [thisEvents; t]; %#ok<AGROW>
        end
    end

    if height(thisEvents) > 0
        thisEvents = sortrows(thisEvents, 'MJD');
    end

    outFile = fullfile(outDir, sprintf('%s_%s.mom', station, compSuffix{c}));
    writeKnownEventsToMom(momFile, thisEvents, compName{c}, outFile);
    fprintf('  -> %s written (%d event(s) within data span)\n', outFile, height(thisEvents));
end

end
