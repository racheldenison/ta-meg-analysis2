function [selectedChannels,D] = meg_runAnalysis(sessionDir)
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

% sessionDir = 'R0817_20181120';

%% session info 

% TA2
exptShortName = 'TA2'; 
exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 
p = meg_params('TA2_Analysis');

% TANoise
% exptShortName = 'TANoise'; % TANoise
% exptDir = '/Users/kantian/Dropbox/Data/TANoise/MEG'; 
% p = meg_params('TANoise_Analysis');

analStr = 'bietfp';
readData = 0; % newly prepare and read data
loadData = 1; % load already prepared data matrix 

selectChannels = 0; 
loadSelectedChannels = 1; 

plotERF = 0; 
plotTF = 0;  
saveTF = 0; % save time frequency mat 

%% setup 

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);
prepDir = sprintf('%s/prep', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));
B = meg_behavior(behav); % update behavior structure

% eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
% eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
% eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

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

data = meg_rejectTrials(data,dataDir); % NaN manually rejected trials

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

behav = meg_behavior(behav); 
conds = behav.responseData_all; % trial x condition 
cond = behav.cuedTarget; % attention condition 

condNames = {'cueCond'}; 

switch exptShortName 
    case 'TA2'
        levelNames = {{'other','cueT1','cueT2','neutral'}}; % levelNames = {{'neutral','cueT1','cueT2','other'}}; 
    case 'TANoise'
        levelNames = {{'other','cueT1','cueT2'}};
end
cond(isnan(cond))=0; 

D = meg_slicer(data, cond, condNames, levelNames); 
% D.cueT2subcueT1 = D.cueT2-D.cueT1; 

field = 'other';
D = rmfield(D,field);

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
% p = [];
p.t = 1:size(data,1);
p.targetWindow = [1000 1400];
classNames = {'V','H'};

for iT = 1:2
    classLabels = B.t1t2Axes(:,iT);
    [A(iT), fH{iT}, figNames{iT}] = meg_plotDecode(data, p, classLabels, classNames);
end

end

 