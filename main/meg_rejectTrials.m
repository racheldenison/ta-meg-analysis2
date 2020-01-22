function data = meg_rejectTrials(sessionDir,data)

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

load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/%s/prep/trials_rejected.mat',sessionDir), 'trials_rejected'); 
data(:,:,trials_rejected) = NaN; 

end

