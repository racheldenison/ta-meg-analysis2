function [A, D, selectedChannels] = meg_runAnalysis(expt, sessionDir, user)

% function [A, D, selectedChannels] = meg_runAnalysis(exptName, sessionDir, [user])
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
    user = 'mcq'; % 'mcq','karen','rachel','karenhd'
end
if ~exist('user','var')
    user = [];
end

%% settings
analStr = 'bietfp';
paramType = 'Analysis';
analType = 'decode'; % 'ERF','TF','decode','channels'
sliceType = 'cue'; % 'all','cue'
channelSelectionType = '20Hz_ebi'; % 'peakprom','classweights','20Hz_ebi'
nChannelsSelected = 5; % number of channels to select from channelsRanked

readData = 0; % need to read data first time
loadData = 1; % reload data matrix
selectChannels = 0;
loadChannels = 1;
flipData = 0; % flip data direction

%%% RD suggestion: switch between these options instead of toggling? so
%%% this function runs only one type of analysis at a time. Suggest
%%% returning an analysis structure A from each plotting function
saveTF = 0; % save time frequency mat

saveAnalysis = 1;
saveFigs = 1;

%% setup
% directories
exptDir = meg_pathToTAMEG(expt, user);
dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);
prepDir = sprintf('%s/prep', dataDir);
figDir = sprintf('%s/figures/%s', dataDir, analStr);
behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);

% file names
fileBase = meg_sessionDirToFileBase(sessionDir, expt);
sqdFile = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file*
dataFile = sprintf('%s/%s_%s_data.mat', matDir, fileBase, analStr);
behavFile = dir(sprintf('%s/*.mat', behavDir));

% eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
% eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase);
% eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase);

% params
p = meg_params(sprintf('%s_%s', expt, paramType));

% behavior
behav = load(sprintf('%s/%s', behavDir, behavFile.name));
B = meg_behavior(behav);

%% make directories

if ~exist(figDir,'dir')
    mkdir(figDir)
end
if ~exist(matDir,'dir')
    mkdir(matDir)
end
if ~exist(prepDir,'dir')
    mkdir(prepDir)
end

%% read data
% read preprocessed sqd, prepare save to prep_data and data matrix
if readData
    [prep_data, data] = meg_getData(sqdFile,p); % time x channels x trials
    save (dataFile, 'data', '-v7.3')
    save (sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
end

% if data already saved, then load data
if loadData
    load(dataFile, 'data')
end

if ~exist('data','var')
    data = [];
end

%% reject trials
data = meg_rejectTrials(data, dataDir); % NaN manually rejected trials

%% channel selection
%%% RD: Reorganizing this section a bit. Let's load sorted channels
%%% (channelsRanked) and then pick some number of top channels from that
%%% list. Channel sorting (making channelsRanked) is now an analysis type
%%% called "channels" in the run analysis section

if loadChannels
    channelsRanked = meg_loadChannels(matDir, channelSelectionType); 
end

selectedChannels = channelsRanked(1:nChannelsSelected); 
disp(sprintf('ch: %d   ',selectedChannels))

%%% RD: Update this part to select nTopChannels. Suggest moving parts of
%%% meg_selectChannels to the new function meg_sortChannels. The current
%%% function meg_selectChannels may become unnecessary, if channel
%%% selection is simply taking the top n channels from the sorted list.

% if selectChannels
%     figPromAvg = sprintf('%spromAvg',figDir);
%     if ~exist(figPromAvg,'dir')
%         mkdir(figPromAvg)
%     end
%     [pkfH, pkfigNames, Pk] = kt_peakstest(sessionDir,data);
%     rd_saveAllFigs(pkfH, pkfigNames, sessionDir, figPromAvg, [])
%     save(sprintf('%s/Pk_avgProm.mat',matDir),'Pk')
%     
%     close all
%     selectedChannels = Pk.passCh';
% else
%     selectedChannels = [];
% end

%% data direction
% flip based on peak prominence direction
if flipData
    data = data.*C.promDir';
end

%% slice data by condition (here cue condition)
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
        
    otherwise
        error('sliceType not recognized')
end

[D, I] = meg_slicer(data, cond, condNames, levelNames);

%% run analysis
switch analType
    case 'channels'
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
        [fH1,figNames1] = meg_plotERF(D,p,selectedChannels);
        
        ERFDir = sprintf('%spromAvg/ERF',figDir);
        if ~exist(ERFDir,'dir')
            mkdir(ERFDir)
        end
        
        rd_saveAllFigs(fH1, figNames1, sessionDir, ERFDir, [])
        
    case 'TF'
        %% TF analyses
        [A,fH2,figNames2] = meg_plotTF(D,p,selectedChannels);
        
        thresholdFigTFDir = sprintf('%spromAvg/TF',figDir);
        if ~exist(thresholdFigTFDir,'dir')
            mkdir(thresholdFigTFDir)
        end
        
        rd_saveAllFigs(fH2, figNames2, sessionDir, thresholdFigTFDir, [])
        if saveTF
            save(sprintf('%s/TF.mat',matDir),'A','-v7.3')
        end
        
    case 'alpha'
        %% eyes closed alpha analyses
        
        [A,fH3,figNames3] = meg_alpha(eyesClosedFileBI,sessionDir);
        rd_saveAllFigs(fH3, figNames3, sessionDir, figDir, [])
        save(sprintf('%s/Alpha.mat',matDir),'A')
        
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

end

