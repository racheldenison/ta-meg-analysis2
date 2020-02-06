function data = meg_rejectTrials(data,dataDir)

% data = MEG_REJECTTRIALS(data)
%
% INPUT
% data
%   data matrix, time x channels x trials
% OUPUT 
% data
%   data matrix, time x channels x trials rejected NaN
%
% Karen Tian
% January 2020

%% NaN manually rejected channels 

load(sprintf('%s/prep/trials_rejected.mat',dataDir), 'trials_rejected'); 
data(:,:,trials_rejected) = NaN; 

end

