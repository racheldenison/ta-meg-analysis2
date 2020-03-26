function [A, selectedChannels, D] = meg_runAnalysis(exptName, sessionDir, user)

% function [A, selectedChannels, D] = meg_runAnalysis(exptName, sessionDir, [user])
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
    exptName = 'TA2'; % TANoise
    sessionDir = 'R0817_20181120';
    user = 'mcq'; % 'mcq','karen'
end
if ~exist('user','var')
    user = [];
end

%% settings 
analStr = 'bietfp';
analType = 'Analysis';

readData = 0; % need to read data first time
loadData = 1; % reload data matrix 
selectChannels = 0; 
loadSelectedChannels = 1; 

%%% RD suggestion: switch between these options instead of toggling? so
%%% this function runs only one type of analysis at a time. Suggest
%%% returning an analysis structure A from each plotting function
plotERF = 0; 
plotTF = 0;  
saveTF = 0; % save time frequency mat 
plotDecode = 0;

%% setup
% i/o
exptDir = meg_pathToTAMEG(exptName, user);
fileBase = meg_sessionDirToFileBase(sessionDir, exptName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);
prepDir = sprintf('%s/prep', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));

% eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
% eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
% eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

% params
p = meg_params(sprintf('%s_%s', exptName, analType));

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
    [prep_data, data] = meg_getData(filename,p); % time x channels x trials
    save (sprintf('%s/%s_%s_data.mat',matDir,fileBase,analStr),'data', '-v7.3')
    save (sprintf('%s/%s_%s_prep_data.mat',prepDir,fileBase,analStr),'prep_data', '-v7.3')
end

% if data already saved, then load data
if loadData
    load(sprintf('%s/%s_%s_data.mat',matDir,fileBase,analStr),'data')
end

%% reject trials
data = meg_rejectTrials(data, matDir); % NaN manually rejected trials

%% threshold channel selector 
% if selectChannels
%     meg_selectChannels(sessionDir,data); 
% elseif loadSelectedChannels 
%     load(sprintf('%s/T.mat',matDir),'T');  
%     selectedChannels = T.passCh; 
% end

if selectChannels
    figPromAvg = sprintf('%spromAvg',figDir);
    if ~exist(figPromAvg,'dir')
        mkdir(figPromAvg)
    end
    [pkfH, pkfigNames, Pk] = kt_peakstest(sessionDir,data);
    rd_saveAllFigs(pkfH, pkfigNames, sessionDir, figPromAvg, [])
    save(sprintf('%s/Pk_avgProm.mat',matDir),'Pk')
    
    close all
    selectedChannels = Pk.passCh';
end

if loadSelectedChannels
    chFile = sprintf('%s/Pk_avgProm.mat',matDir); 
    load(chFile)
    selectedChannels = Pk.passCh';
    channelDir = Pk.promDir(selectedChannels); % positive or negative peak 
end

%% data direction 
% flip based on peak prominence direction 
data = data.*Pk.promDir'; 

%% peak channel selector

% try peak channels 
% [fHpeak,figNamesPk] = meg_selectChannels(sessionDir, data); 
% figDirPeak = sprintf('%speak',figDir);
% if ~exist(figDirPeak,'dir')
%     mkdir(figDirPeak)
% end
% rd_saveAllFigs(fHpeak, figNamesPk, sessionDir, figDirPeak, [])
% close all

%% try alpha channel selector 
% C = meg_selectChannels(sessionDir); 

% selectedChannels = C.chSortAlpha(1:5)'; 

% figDirAlpha = sprintf('%sAlphaTop5',figDir);
% if ~exist(figDirAlpha,'dir')
%     mkdir(figDirAlpha)
% end

%% slice data by condition (here cue condition) 
cond = B.cuedTarget; % attention condition 
condNames = {'cueCond'}; 

switch exptName 
    case 'TA2'
        levelNames = {{'cueT1','cueT2','neutral'}}; % levelNames = {{'neutral','cueT1','cueT2','other'}}; 
    case 'TANoise'
        levelNames = {{'cueT1','cueT2'}};
end

[D, I] = meg_slicer(data, cond, condNames, levelNames); 
% D.cueT2subcueT1 = D.cueT2-D.cueT1; 

%% ERF analyses

if plotERF
    [fH1,figNames1] = meg_plotERF(D,p,selectedChannels);
    
    ERFDir = sprintf('%spromAvg/ERF',figDir);
    if ~exist(ERFDir,'dir')
        mkdir(ERFDir)
    end
    
    rd_saveAllFigs(fH1, figNames1, sessionDir, ERFDir, [])
    close all
end

%% TF analyses

if plotTF
    [TF,fH2,figNames2] = meg_plotTF(D,p,selectedChannels);
    
    thresholdFigTFDir = sprintf('%spromAvg/TF',figDir);
    if ~exist(thresholdFigTFDir,'dir')
        mkdir(thresholdFigTFDir)
    end
    
    rd_saveAllFigs(fH2, figNames2, sessionDir, thresholdFigTFDir, [])
    if saveTF
        save(sprintf('%s/TF.mat',matDir),'TF','-v7.3')  
    end
    close all
end

%% eyes closed alpha analyses 

% [Alpha,fH3,figNames3] = meg_alpha(eyesClosedFileBI,sessionDir);
% rd_saveAllFigs(fH3, figNames3, sessionDir, figDir, [])
% save(sprintf('%s/Alpha.mat',matDir),'Alpha')

%% Decoding orientation
if plotDecode
    % p = [];
    p.t = 1:size(data,1);
    p.targetWindow = [1000 1400];
    classNames = {'V','H'};
    
    for iT = 1:2
        classLabels = B.t1t2Axes(:,iT);
        [A(iT), fH{iT}, figNames{iT}] = meg_plotDecode(D, I, p, classLabels, classNames);
    end
end

end

 