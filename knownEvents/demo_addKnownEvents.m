%% demo_addKnownEvents.m
% Demonstrates the known-offset import pipeline on the sample files
% provided alongside this package:
%   HKSC_0.mom, HKSC_1.mom, HKSC_2.mom  - target N/E/U mom files
%   steps.txt                            - NGL steps catalogue (HKSC has
%                                           entries -> Hong Kong stations
%                                           ARE tracked by NGL)
%   BJFS.txt                              - CEA header-annotated coordinate
%                                           file (Mainland China stations,
%                                           NOT in NGL, use this path
%                                           instead of/alongside StepsFile)

dataDir = fileparts(mfilename('fullpath'));

% --- Case 1: station available in NGL (e.g. HKSC) -----------------------
addKnownEventsToMom('HKSC', dataDir, fullfile(dataDir,'out'), ...
    'StepsFile', fullfile(dataDir,'steps.txt'));

% --- Case 2: Mainland China station, NGL has nothing, use CEA header ----
% (BJFS mom files are not part of this sample set; this call only shows
%  the header parser output. In real use, point momDir at BJFS_0/1/2.mom.)
T_bjfs = parseCEAHeader(fullfile(dataDir,'BJFS.txt'));
disp(T_bjfs);

% --- Case 3: a China station with BOTH sources available -----------------
% addKnownEventsToMom('XXXX', momDir, outDir, ...
%     'StepsFile', 'steps.txt', 'CEAFile', 'XXXX.txt');
