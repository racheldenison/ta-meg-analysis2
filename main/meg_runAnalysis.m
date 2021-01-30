function [A, D, selectedChannels, I, B] = meg_runAnalysis(expt, sessionDir, user)

% function [A, D, selectedChannels, I, B] = meg_runAnalysis(exptName, sessionDir, [user])
%
% Loads behavioral data, loads MEG data, rejects trials, selects
% channels, slices MEG data by condition, does analysis (ERF, TF, ITPC,
% decoding). 
% 
% INPUTS
% exptName
%   string, 'TA2' or 'TANoise'
%
% sessionDir
%   string, e.g., 'R0817_20181120'
%
% user (optional)
%   string giving the user name to determine path to data. defaults to
%   'mcq' = get the data from the mcq server
%
% OUTPUTS 
% A 
%   Analysis structure 
% D
%   Data structure sliced by condition
%   Each field is a condition combination containing a data matrix (time x 
%   channels x trials) for the trials in that group.
% I
%   trial indices for each condition combination
% B
%   Structure of behaviorial data (including cued target, cue validity,
%   response target, target orientation, accuracy, reaction time, target axes
%   (horizontal/vertical), ITI jitter)
%
% Karen Tian
% January 2020

%% inputs
if nargin == 0
    expt = 'TA2'; % 'TANoise'
    sessionDir = 'R1187_20181119'; % 'R1187_20180105'
    user = 'karen'; % 'mcq','karen','rachel','karenhd'
end
if ~exist('user','var')
    user = [];
end

%% Settings
% analysis options
analStr = 'ebi'; % 'bietfp'
preprocStr = 'ebi'; 

paramType = 'ITPCsession8'; % for meg_params 'Preproc', 'Analysis', 'ITPC', 'ITPCsession8" (to correct for the cutoff run)
analType = 'none'; % 'ITPC'; % 'none','readdata','sortchannels','ERF','TF','decode', 'TFwholeTrial'
avgTrial = 0; 
sliceType = 'cue'; % 'all','cue','cueAcc','ITIjitter','ITICue'

% data
getData = 'fromSqd'; % 'fromSqd' (read data from sqd), 'fromMat' (load data from prepared .mat) 
saveData = 0; % save .mat 

% channels
loadChannels = 1;
channelSelectionType = '20Hz_ebi'; % '20Hz_ebi' (rank by 20Hz power), 'peakprom' (rank by T1 and T2 ERP peak prominence), 'classweights'

selectChannels = 1;
nChannelsSelected = 5; % number of channels to select from channelsRanked

% saving
saveAnalysis = 0;
saveFigs = 1;

%% Setup
% directories
exptDir = meg_pathToTAMEG(expt, user);
dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc_%s', dataDir, preprocStr);
prepDir = sprintf('%s/prep', dataDir);
figDir = sprintf('%s/figures/%s', dataDir, analStr);
behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
stimDir = sprintf('%s/Behavior/%s/stimuli', exptDir(1:end-4), sessionDir);

% file names
fileBase = meg_sessionDirToFileBase(sessionDir, expt);
sqdFile = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file*
dataFile = sprintf('%s/%s_%s_data.mat', matDir, fileBase, analStr);
behavFile = dir(sprintf('%s/*.mat', behavDir));

% params
p = meg_params(sprintf('%s_%s', expt, paramType));

% behavior
behav = load(sprintf('%s/%s', behavDir, behavFile.name));
B = meg_behavior(behav);

% stimuli 
stimFiles = dir(sprintf('%s/*.*',stimDir));
stimFiles(2) = []; stimFiles(1) = []; 
for i = 1:numel(stimFiles)
    stimFileNames{i} = stimFiles(i).name; 
end
sortedFiles = natsortfiles(stimFileNames); 
jitSeq = []; 
for i = 1:numel(stimFiles) % 1:numel(stimFiles) % for one subject, 2:13
    % stimFileName{i} = fullfile(stimDir,stimFiles(i).name);
    % stimFileName = natsortfiles(stimFileName); 
    stim = load(sprintf('%s/%s',stimDir,sortedFiles{i})); 
    jitSeq = [jitSeq stim.stimulus.jitSeq]; 
end
B.jitSeq = jitSeq'; % jitters in s 
jitters = unique(jitSeq);
condVal = 1; % replace s with int 
for i = 1:numel(jitters)
    B.jitSeqCond(B.jitSeq==jitters(i)) = condVal; 
    condVal = condVal + 1; 
end
B.jitSeqCond = B.jitSeqCond'; 
B.jitSeqCond(isnan(B.targetPresent)) = NaN;

%% Make directories
if ~exist(figDir,'dir')
    mkdir(figDir)
end
if ~exist(matDir,'dir')
    mkdir(matDir)
end
if ~exist(prepDir,'dir')
    mkdir(prepDir)
end

%% Get data (time x ch x trials) 
switch getData        
    case 'fromSqd' % Read preprocessed sqd using fieldtrip 
        [~, data] = meg_getData(sqdFile,p); % time x channels x trials
        if saveData
            save(dataFile, 'data', '-v7.3')
            save(sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
        end
    case 'fromMat'
        load(dataFile, 'data') % Load time x ch x trial mat file 
    otherwise
        data = [];
end
fprintf('data read %s \n', getData)

%% Reject trials
[data, nTrialsRejected] = meg_rejectTrials(data, dataDir); % NaN manually rejected trials
fprintf('%d trials rejected \n', nTrialsRejected)

%% flip data based on amplitude at precue 
% switch analType
%     case 'ITPC'
%         load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/channels/TANoise_20Hz_chDir.mat','chDir')
%         chDir = chDir(sessionIdx,:); % resave w sessionDir info so don't need sessionIdx input 
%         
%         for iCh = 1:numel(p.megChannels)
%             vals = [];
%             flip = 1;
%             if chDir(iCh)==-1
%                 flip = -1;
%             end
%             vals = data(:,iCh,:);
%             vals = vals.*flip;
%             data(:,iCh,:) = vals;
%         end
%         
%         disp('flipped data based on polarity')
%         
%     otherwise
%         disp('unflipped data')
%         
% end

%% Channel selection
if loadChannels
    if strcmp(channelSelectionType,'20Hz_ebi') && strcmp(expt,'TA2')
        error('For TA2, change channelSelectionType to peakprom') 
    end
    channelsRanked = meg_loadChannels(matDir, channelSelectionType); 
end

if selectChannels
    selectedChannels = channelsRanked(1:nChannelsSelected);
    disp(sprintf('ch: %d   ',selectedChannels))
else
    selectedChannels = [];
end

fprintf('channel selection: %s \n', channelSelectionType)

%% Slice data by condition 
switch sliceType
    case 'all'
        cond = B.targetPresent;
        condNames = {'trial'};
        levelNames = {{'all'}};
        
    case 'cue'
        cond = B.cuedTarget;
        condNames = {'cueCond'};
        switch expt
            case 'TA2'
                levelNames = {{'cueT1','cueT2','neutral'}};
            case 'TANoise'
                levelNames = {{'cueT1','cueT2'}};
        end
    
    case 'cueAcc' % cue, accruacy 
        cond = [B.cuedTarget B.acc]; 
        condNames = {'cueCond','acc'}; 
        switch expt 
            case 'TA2'
                levelNames = {{'cueT1','cueT2','cueN'},{'incorrect','correct'}}; 
           case 'TANoise'
                levelNames = {{'cueT1','cueT2'},{'incorrect','correct'}}; 
        end 
               
    case 'ITIjitter'
        if strcmp(expt,'TA2')
            error('Slice type ITI jitter for TANoise only. ITI jitter determines onset of 20Hz flicker')
        end
        cond = B.jitSeqCond; 
        condNames = {'ITIjitter'}; 
        levelNames = {sprintfc('ITI%d',500:200:1500)};
        
    case 'ITICue'
        if strcmp(expt,'TA2')
            error('Slice type ITI jitter for TANoise only. ITI jitter determines onset of 20Hz flicker')
        end
        cond = [B.acc B.jitSeqCond];
        condNames = {'cueCond','ITIjitter'}; 
        levelNames = {{'cueT1','cueT2'},sprintfc('ITI%d',500:200:1500)};
        
    otherwise 
        error('sliceType not recognized')
end

[D, I] = meg_slicer(data, cond, condNames, levelNames);
fprintf('sliced data by: %s \n', sliceType)

%% Run analysis
switch analType
    case 'none'       
        A = []; % just return sliced data or selected channels
        
    case 'readdata'
        %% Read preprocessed sqd
        [prep_data, data] = meg_getData(sqdFile,p); % time x channels x trials
        save (dataFile, 'data', '-v7.3')
        save (sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
        
    case 'sortchannels'
        %% Sort channels
        %%% RD: Update meg_sortChannels
        sortType = 'classweights';
        
        [C, fH, figNames] = meg_sortChannels(expt, sessionDir, user, sortType);
        A = C;
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, [], figDir)
        end
        if saveAnalysis
             save(sprintf('%s/channels_%s.mat', matDir, sortType), 'C')
        end
        
    case 'ERF'
        %% ERF analyses
        [A,fH,figNames] = meg_plotERF(D,p,selectedChannels);
        
        ERFDir = sprintf('%spromAvg/ERF',figDir);
        if ~exist(ERFDir,'dir')
            mkdir(ERFDir)
        end
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, ERFDir, [])
        end
        
    case 'TF'
        %% TF analyses
        if avgTrial
            fieldNames = fieldnames(D); Davg = [];
            for i = 1:numel(fieldNames)
                Davg.(fieldNames{i})(:,:,1) = nanmean(D.(fieldNames{i}),3);
            end
            D = []; D = Davg;
        end
        [A,fH,figNames] = meg_plotTF(D,p,selectedChannels);
        
        thresholdFigTFDir = sprintf('%s/TF_singleTrial_condCue',figDir);
        if ~exist(thresholdFigTFDir,'dir')
            mkdir(thresholdFigTFDir)
        end
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, thresholdFigTFDir, [])
        end
        if saveAnalysis
            % save(sprintf('%s/TF.mat',matDir),'A','-v7.3')
            save(sprintf('%s/TFspectrogram_5Ch_avgTrial.mat',matDir),'A','-v7.3')
        end
        
    case 'TFwholeTrial'
        %% TF analyses
        selectedChannels = channelsRanked(1:5);  

        if avgTrial 
          fieldNames = fieldnames(D); Davg = [];
            for i = 1:numel(fieldNames)
                Davg.(fieldNames{i})(:,:,1) = nanmean(D.(fieldNames{i}),3);
            end
            D = []; D = Davg;
        end
        [A,fH,figNames] = meg_plotFFT_wholeTrial(D,p,selectedChannels); 
        
        % thresholdFigTFDir = sprintf('%s/TF_wholeTrial_5Ch',figDir);
        % thresholdFigTFDir = sprintf('%s/TF_prePreCue_5Ch',figDir);
        thresholdFigTFDir = sprintf('%s/TF_5Ch_averageTrial',figDir);
        if ~exist(thresholdFigTFDir,'dir')
            mkdir(thresholdFigTFDir)
        end
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, thresholdFigTFDir, [])
        end
        if saveAnalysis
            % save(sprintf('%s/TF_wholeTrial.mat',matDir),'A','-v7.3')
            % save(sprintf('%s/TF_prePreCue.mat',matDir),'A','-v7.3')
            save(sprintf('%s/TF_preTarget.mat',matDir),'A','-v7.3')
        end
        
    case 'ITPC' 
        %% time freq and ITPC phase angle 
        % selectedChannels = p.megChannels; % all channels 
        % selectedChannels = channelsRanked(1:5); % top 5 channels 
        if avgTrial
            % average data across trials
            fieldNames = fieldnames(D); Davg = [];
            for i = 1:numel(fieldNames)
                Davg.(fieldNames{i})(:,:,1) = nanmean(D.(fieldNames{i}),3);
            end
            D = []; D = Davg;
            [A,fH,figNames] = meg_plotITPCavg(D,sessionDir,p,selectedChannels,p.ssvefFreq); 
        else % single trial 
            % selectedFreq = p.ssvefFreq; 
            selectedFreq = 1:50; 
            [A,fH,figNames] = meg_plotITPC(D,sessionDir,p,selectedChannels,selectedFreq); 
        end
        
        itpcFigDir = sprintf('%s/ITPCspectrogram_singleTrial_ITICue',figDir); 
        if ~exist(itpcFigDir,'dir')
            mkdir(itpcFigDir)
        end
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, itpcFigDir, [])
        end
        if saveAnalysis
            save(sprintf('%s/ITPCspectrogram.mat',matDir),'A','-v7.3')
        end
        
    case 'alpha'
        %% Eyes closed alpha analyses
        % RD: adjust inputs, rename to meg_plotAlpha? revisit 
        [A,fH,figNames] = meg_alpha(eyesClosedFileBI,sessionDir);
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])
        end
        if saveAnalysis
            save(sprintf('%s/Alpha.mat',matDir),'A')
        end
        
    case 'decode'       
        %% Decoding orientation
        analysisName = 'classAcc';
        targetNames = {'T1','T2'};
        classNames = {'V','H'};
        twin = [-50 550];
        condNames = fields(D);
        
%         allA = [];
        nTopCh = [157 100 50 20 10];
        i = 3;
%         for i = 1:numel(nTopCh)
            
            selectedChannels = channelsRanked(1:nTopCh(i));
            
%             A = struct([]);
            for iT = 1:2
                targetTime = p.eventTimes(strcmp(p.eventNames,targetNames{iT}));
                p.targetWindow = targetTime + twin;
                classLabels = B.t1t2Axes(:,iT);
                [A(iT), fH{iT}, figNames{iT}] = meg_plotDecode(D, I, p, classLabels, classNames, selectedChannels);
                
                if saveFigs
                    for iF = 1:numel(fH{iT})
                        figNamesT{iF} = sprintf('%s_%s', figNames{iT}{iF}, targetNames{iT});
                    end
                    rd_saveAllFigs(fH{iT}, figNamesT, [], figDir)
                end
            end
            
%             allA{i} = A;
%         end
        
        if saveAnalysis
            decodeAnalStr = A(1).(condNames{1}).decodingOps.analStr;
            save(sprintf('%s/%s_%s_%sSlice_%s.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'A', 'targetNames')
%             save(sprintf('%s/%s_%s_%sSlice_%s_varyNChannels10-157_sortedBy20Hz.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'allA', 'targetNames', 'nTopCh')
            save(sprintf('%s/%s_%s_%sSlice_%s_nCh%d.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr, nTopCh(i)), 'A', 'targetNames')
        end
        
    otherwise
        error('analType not recognized')
end
fprintf('analysis: %s \n', analType)

end

