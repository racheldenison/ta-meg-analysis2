%% load group 

%% setup
user = 'karen'; % 'mcq','karen','rachel'
expt = 'TANoise'; % 'TANoise'
sessionIdx = 1:20; 
p = meg_params('TANoise_ITPCsession8'); 

%% get session names
allSessions = meg_sessions(expt);
sessionNames = allSessions(sessionIdx);

%% run analysis 
groupA = []; groupD = []; groupI = []; groupB = []; 

%%
for i = 2:numel(sessionNames)
    sessionDir = sessionNames{i}; 
    
    %% setup
    analStr = 'ebi'; % 'bietfp'
    preprocStr = 'ebi';
    
    paramType = 'ITPCsession8'; % for meg_params 'Preproc', 'Analysis', 'ITPC', 'ITPCsession8" (to correct for the cutoff run)
    analType = 'ITPC'; % 'none','readdata','sortchannels','ERF','TF','decode'
    avgTrial = 0;
    sliceType = 'all'; % 'all','cue','cueAcc'

    exptDir = meg_pathToTAMEG(expt, user);
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    preprocDir = sprintf('%s/preproc_%s', dataDir, preprocStr);
    prepDir = sprintf('%s/prep', dataDir);
    figDir = sprintf('%s/figures/%s', dataDir, analStr);
    behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
    % file names
    fileBase = meg_sessionDirToFileBase(sessionDir, expt);
    sqdFile = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file*
    dataFile = sprintf('%s/%s_%s_data.mat', matDir, fileBase, analStr);
    behavFile = dir(sprintf('%s/*.mat', behavDir));
    % behavior
    behav = load(sprintf('%s/%s', behavDir, behavFile.name));
    B = meg_behavior(behav);
    
    %% fileOI
    load(sprintf('%s/ITPC.mat',matDir)); % ITPC mat A with 157 channels 
    
    groupA(i).tfAmps = A.all.tfAmps;
    groupA(i).tfPpws = A.all.tfPows; 
    groupA(i).ITPC = A.all.ITPC;
end

