function [figDir,dateStr,style,colors,p] = meg_manuscriptParams
% function [opt,figDir,style,colors] = meg_manuscriptParams
% Setup figure directory
% Load ITPC params and figure style and colors 

%% 
[style, colors] = meg_manuscriptStyle; % load figure styling and colors 
p = meg_params('TANoise_ITPCsession8'); % load ITPC params 

%% Figure directory
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
user = 'kantian'; % karen 
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/manuscriptFigures/figs',user);
if ~exist(figDir, 'dir')
    mkdir(figDir)
end
