function meg_runAnalysis(sessionDir)
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

%% session info 

sessionDir = 'R0817_20181120';
exptShortName = 'TA2'; % TANoise
analStr = 'bietfp';
% exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 
exptDir = meg_pathToTA2('MEG');

p = meg_params('TA2_Analysis');
readData = 0; % need to read data first time
loadData = 1; % reload data matrix 

%% setup 

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));
B = meg_behavior(behav); % update behavior structure

eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

%% make directories

if ~exist(figDir,'dir')
    mkdir(figDir)
end
if ~exist(matDir,'dir')
    mkdir(matDir)
end

%% read data

if readData
    data = meg_getData(filename,p); % time x channels x trials
    % save data
    save(sprintf('%s/data.mat',matDir),'data')
end

% if data already saved, then load data
if loadData
    load(sprintf('%s/data.mat',matDir),'data')
end

%% reject trials

data = meg_rejectTrials(data,dataDir); % NaN rejected trials

%% peak channel selector

% try peak channels 
[fHpeak,figNamesPk] = meg_selectChannels(sessionDir, data); 
figDirPeak = sprintf('%speak',figDir);
if ~exist(figDirPeak,'dir')
    mkdir(figDirPeak)
end
rd_saveAllFigs(fHpeak, figNamesPk, sessionDir, figDirPeak, [])
close all

%% try alpha channel selector 
% C = meg_selectChannels(sessionDir); 

% selectedChannels = C.chSortAlpha(1:5)'; 
% 
% figDirAlpha = sprintf('%sAlphaTop5',figDir);
% if ~exist(figDirAlpha,'dir')
%     mkdir(figDirAlpha)
% end

%% slice data by condition (here cue condition) 

% behav = meg_behavior(behav); 
% conds = behav.responseData_all; % trial x condition 
% cond = behav.cuedTarget; % attention condition 
% 
% condNames = {'cueCond'}; 
% levelNames = {{'other','cueT1','cueT2'}}; 
% 
% D = meg_slicer(data, cond, condNames, levelNames); 
% D.cueT2subcueT1 = D.cueT2-D.cueT1; 

%% ERF analyses

% [fH1,figNames1] = meg_plotERF(D,selectedChannels); 
% rd_saveAllFigs(fH1, figNames1, sessionDir, figDirAlpha, [])
% close all

%% TF analyses

% [A,fH2,figNames2] = meg_plotTF(D,selectedChannels); 
% rd_saveAllFigs(fH2, figNames2, sessionDir, figDirAlpha, [])
% close all

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

 