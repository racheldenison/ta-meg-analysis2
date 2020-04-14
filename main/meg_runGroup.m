% meg_runGroup 

%% setup
user = 'mcq'; % 'mcq','karen','rachel'
expt = 'TANoise';
sessionIdx = 1:20; 

%% get session names
allSessions = meg_sessions(expt);
sessionNames = allSessions(sessionIdx);

%% run analysis 
groupA = [];
for i=1:numel(sessionNames) 
    sessionDir = sessionNames{i}; 
    disp(sessionDir)

    [A, D, selectedChannels] = meg_runAnalysis(expt, sessionDir, user); 
    groupA{i} = A;
    pause(1)
    close all
end