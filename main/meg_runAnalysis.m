% meg_runAnalysis
% Karen Tian 
% January 2020 

%% session info 

sessionDir = 'R1507_20190627';
exptShortName = 'TA2';
analStr = 'bietfp';
exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 

p = meg_params('TA2');

%% setup 

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
% savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBa se, analStr);
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));

%% read data

data = meg_getData(filename,p); 

%% reject trials

load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/%s/prep/trials_rejected.mat',sessionDir)); 
data(:,:,trials_rejected) = NaN; 

%% analyses 


% feed to slicer 
% channel selection
    % peak, snr
    % plot on topo L/R 
    % 1,3,5,7

fH = meg_plotERF(); 
[A,fH2] = meg_plotTF(); 





