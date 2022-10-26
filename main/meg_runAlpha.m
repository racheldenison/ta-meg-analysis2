% meg_runAlpha
%% input 
expt = 'TA2'; % 'TANoise' 'TA2'
sessionNames = meg_sessions(expt);
user = 'karen';
saveFigs = 1; 
saveAnalysis = 1; 

%%
for i = 1:numel(sessionNames)
    sessionDir = sessionNames{i}; % 'R0817_20181120';

    %% setup
    % directories
    exptDir = meg_pathToTAMEG(expt, user);
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    prepDir = sprintf('%s/prep', dataDir);
    figDir = sprintf('%s/figures/eyesClosed', dataDir);
    
    if ~exist(figDir,'dir')
        mkdir(figDir)
    end
    
    % file names
    fileBase = meg_sessionDirToFileBase(sessionDir, expt);
    eyesShortName = strrep(fileBase,expt,'EyesClosed');
    eyesClosedFile = sprintf('%s/%s.sqd',dataDir,eyesShortName);
    
    %% alpha analysis
    [A,fH,figNames] = meg_alpha(eyesClosedFile,sessionDir);
    
    if saveFigs
        rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])
    end
    if saveAnalysis
        save(sprintf('%s/Alpha.mat',matDir),'A')
    end
    close all
end
