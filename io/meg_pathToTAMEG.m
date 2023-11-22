function path = meg_pathToTAMEG(exptName, user)
% function path = meg_pathToTAMEG(exptName, user)

if nargin < 1
    error('Must provide experiment name [TANoise or TA2]')
end
if nargin < 2 || isempty(user)
    user = 'mcq';
end

switch user
    case 'mcq'
        root = '/Volumes/purplab2/EXPERIMENTS/1_Current_Experiments/Rachel/TA_MEG/MEG';
    case 'karen'
        root = '/Users/kantian/Dropbox/Data';
    case 'karenhd'
        root = 'Volumes/kantian/Dropbox/Data'; 
    case 'rachel'
        root = '/Volumes/Rachel NYU/Backup_mcq_20200321/TA_MEG/MEG';
    case 'scc'
        root = '/projectnb/rdenlab';
    otherwise
        error('user not recognized')
end

path = sprintf('%s/%s/MEG', root, exptName);
