function [A, D, selectedChannels, I, B] = meg_runAnalysis(expt, sessionDir, user,sessionIdx)

% function [A, D, selectedChannels, I, B] = meg_runAnalysis(exptName, sessionDir, [user])
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
% %%% RD: Update documentation about what this function can do
% setup
% read data
% reject trials
% select channels
% slicer
% analysis: ERF
% analysis: TF
% analysis: alpha eyes closed
%
% Karen Tian
% January 2020


%% inputs
if nargin == 0
    expt = 'TA2'; % TANoise
    sessionDir = 'R1187_20181119';
%     expt = 'TANoise';
%     sessionDir = 'R1187_20180105';
    user = 'karen'; % 'mcq','karen','rachel','karenhd'
end
if ~exist('user','var')
    user = [];
end

%% Settings
% analysis options
analStr = 'ebi'; % 'bietfp'
preprocStr = 'ebi'; 

paramType = 'ITPC'; % for meg_params 'Preproc', 'Analysis', 'ITPC', 'ITPCsession8" (to correct for the cutoff run)
analType = 'ITPC'; % 'none','readdata','sortchannels','ERF','TF','decode'
sliceType = 'all'; % 'all','cue','cueAcc'

% data
% readData = 0; % read data from sqd into matrix .mat 
% loadData = 1; % read data from prepared .mat 
getData = 'fromSqd'; % 'fromSqd' (read data from sqd), 'fromMat' (load data from prepared .mat) 
saveData = 0; % save .mat 

% channels
loadChannels = 1;
channelSelectionType = '20Hz_ebi'; % '20Hz_ebi'; % 'peakprom','classweights','20Hz_ebi'

selectChannels = 0;
nChannelsSelected = 5; % number of channels to select from channelsRanked

% saving
saveAnalysis = 1;
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
    case 'fromSqd' % Read preprocessed sqd
        [~, data] = meg_getData(sqdFile,p); % time x channels x trials
        if saveData
            save(dataFile, 'data', '-v7.3')
            save(sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
        end
    case 'fromMat'
        load(dataFile, 'data') % time x ch x trial 
    otherwise
        data = [];
end
fprintf('data read %s \n', getData)

%% Reject trials
[data, nTrialsRejected] = meg_rejectTrials(data, dataDir); % NaN manually rejected trials
fprintf('%d trials rejected \n', nTrialsRejected)

%% flip data based on polarity

switch analType
    case 'ITPC'
        load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/channels/TANoise_20Hz_chDir.mat','chDir')
        chDir = chDir(sessionIdx,:); % resave w sessionDir info so don't need sessionIdx input 
        
        for iCh = 1:numel(p.megChannels)
            vals = [];
            flip = 1;
            if chDir(iCh)==-1
                flip = -1;
            end
            vals = data(:,iCh,:);
            vals = vals.*flip;
            data(:,iCh,:) = vals;
        end
        
        disp('flipped data based on polarity')
        
    otherwise
        disp('unflipped data')
        
end

%% Channel selection
if loadChannels
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
    
    case 'cueAcc'
        cond = [B.cuedTarget B.acc]; 
        condNames = {'cueCond','acc'}; 
        levelNames = {{'cueT1','cueT2','cueN'},{'incorrect','correct'}}; 
        % cond(acc==1) = cond(acc==1)+10; % correct      
               
    otherwise 
        error('sliceType not recognized')
end

[D, I] = meg_slicer(data, cond, condNames, levelNames);
fprintf('sliced data by: %s \n', sliceType)

%% Run analysis
switch analType
    case 'none'
        % just return sliced data or selected channels
        A = [];
        
%     case 'readdata'
%         %% Read preprocessed sqd
%         [prep_data, data] = meg_getData(sqdFile,p); % time x channels x trials
%         save (dataFile, 'data', '-v7.3')
%         save (sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
        
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
        [A,fH,figNames] = meg_plotTF(D,p,selectedChannels);
        
        thresholdFigTFDir = sprintf('%spromAvg/TF',figDir);
        if ~exist(thresholdFigTFDir,'dir')
            mkdir(thresholdFigTFDir)
        end
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, thresholdFigTFDir, [])
        end
        if saveAnalysis
            save(sprintf('%s/TF.mat',matDir),'A','-v7.3')
        end
        
    case 'ITPC' % phase angle 
        
        % average data across trials 
        fieldNames = fieldnames(D); Davg = []; 
        for i = 1:numel(fieldNames)
            Davg.(fieldNames{i})(:,:,1) = nanmean(D.(fieldNames{i}),3);
        end
        D = []; D = Davg; 
        
        selectedChannels = channelsRanked(1:5); 
        [A,fH,figNames] = meg_plotITPCavg(D,sessionDir,p,selectedChannels,p.ssvefFreq); 
        
        itpcFigDir = sprintf('%s/ITPCAvgTrials',figDir); 
        if ~exist(itpcFigDir,'dir')
            mkdir(itpcFigDir)
        end
        
        if saveFigs
            rd_saveAllFigs(fH, figNames, sessionDir, itpcFigDir, [])
        end
        if saveAnalysis
            save(sprintf('%s/ITPCavg.mat',matDir),'A','-v7.3')
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
%             save(sprintf('%s/%s_%s_%sSlice_%s.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'A', 'targetNames')
%             save(sprintf('%s/%s_%s_%sSlice_%s_varyNChannels10-157_sortedBy20Hz.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'allA', 'targetNames', 'nTopCh')
            save(sprintf('%s/%s_%s_%sSlice_%s_nCh%d.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr, nTopCh(i)), 'A', 'targetNames')
        end
        
    otherwise
        error('analType not recognized')
end

fprintf('analysis: %s \n', analType)

end

