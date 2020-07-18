function [data,nTrialsRejected] = meg_rejectTrials(data, dataDir)

% data = MEG_REJECTTRIALS(data)
%
% INPUT
% data
%   data matrix, time x channels x trials
% OUPUT 
% data
%   data matrix, time x channels x trials rejected NaN
% nTrialsRejected
%   number of trials rejected 
%
% Karen Tian
% January 2020

%% check inputs
if isempty(data)
    return
end

%% NaN manually rejected channels 
try
    load(sprintf('%s/mat/trials_rejected.mat', dataDir), 'trials_rejected'); % file paths of trials_rejected sometimes saved differently 
catch
    load(sprintf('%s/prep/trials_rejected.mat', dataDir), 'trials_rejected'); 
end
nTrialsRejected = size(trials_rejected,1); 
data(:,:,trials_rejected) = NaN; 

end

