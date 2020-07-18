% meg_runGroup 
% group meg_runAnalysis 

%% setup
user = 'karen'; % 'mcq','karen','rachel'
expt = 'TANoise'; % 'TANoise'
sessionIdx = 1:20; 

%% get session names
allSessions = meg_sessions(expt);
sessionNames = allSessions(sessionIdx);

%% run analysis 
groupA = []; groupD = []; groupI = []; groupB = []; 

for i= 2:7% 1:numel(sessionNames) % skipped 8, cant read sqd? m
    sessionDir = sessionNames{i}; 
    disp(sessionDir)
    sessionIdx = i; 

    [A, D, selectedChannels, I, B] = meg_runAnalysis(expt,sessionDir,user,sessionIdx); 
    % [A] = loadGroupAnalysis(expt,sessionDir,user); 
    groupA{i} = A;
    % groupD{i} = D; 
    % groupB{i} = B; 
    close all
end

%% average across trials
fields = fieldnames(groupD{1}); 
for i=1:numel(sessionNames)
    for iF=1:numel(fields) % incorrect trials only 
    vals = nanmean(groupD{i}.(fields{iF}),3); 
    groupDavg.(fields{iF})(:,:,i) = vals; 
    end
end

%% compile group (keep trials)
fields = fieldnames(groupD{1});
for i=1:numel(sessionNames) 
    for iF=1:numel(fields)
        disp(sessionNames(i)) 
        vals = []; 
        vals = groupD{i}.(fields{iF});
        groupDavg.(fields{iF})(:,:,:,i) = vals;
    end
end






