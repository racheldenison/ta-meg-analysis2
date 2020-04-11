function data = meg_rejectTrials(data, dataDir)

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

%% check inputs
if isempty(data)
    return
end

%% NaN manually rejected channels 
try
    load(sprintf('%s/mat/trials_rejected.mat', dataDir), 'trials_rejected'); 
catch
    load(sprintf('%s/prep/trials_rejected.mat', dataDir), 'trials_rejected'); 
end
data(:,:,trials_rejected) = NaN; 

end

