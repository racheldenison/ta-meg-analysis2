function C = meg_selectChannels(sessionDir)

% MEG_SELECTCHANNELS(sessionDir)
%
% INPUTS
% sessionDir 
% 
% OUPUT
%   C
%       structure of top channels by analysis type 
%
% Karen Tian
% January 2020

% snr, peak, topo, 

%% setup

% sessionDir = 'R0817_20190625';

exptShortName = 'TA2';
analStr = 'bietfp';
exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 

p = meg_params('TA2');

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));

eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

%% select based off alpha 

load(sprintf('%s/Alpha.mat',matDir)); 

power = Alpha.TFRhann.powspctrm; 
powerAvg = nanmean(power,3); % average power across time 
alpha = nanmean(powerAvg(:,Alpha.tfAlpha-2:Alpha.tfAlpha+2),2); % individual alpha ± 2Hz, avg across time
[C.powSortAlpha,C.chSortAlpha] = sort(alpha,'descend'); 

figure
set(gcf,'Position',[100 100 300 500])
imagesc(alpha); % check
set(gca,'xtick',[],'ytick',0:10:160)
ylabel('channel')
title(sprintf('%s eyes closed alpha %d ± 2Hz \n top 5 channels [%d %d %d %d %d]',und2space(sessionDir),Alpha.tfAlpha,C.chSortAlpha(1:5)))

%% save figs

fH = sort(findobj('Type','figure'));
figNames = {'TopChannels_Alpha'};
rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])

%% save top channels 

save(sprintf('%s/C.mat',matDir),'C')

end



