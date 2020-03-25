function data = meg_rejectTrials(data,trDir)

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

load(sprintf('%s/trials_rejected.mat',trDir), 'trials_rejected'); 
data(:,:,trials_rejected) = NaN; 

end

