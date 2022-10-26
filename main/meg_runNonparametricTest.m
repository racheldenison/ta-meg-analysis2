function [A] = meg_runNonparametricTest(data1,data2,nPerm,stat)

%% 
nPerm = 1000; 
stat = 'sum'; 

data1 = A.cueT1.normSubjectFlipped;
data2 = A.cueT2.normSubjectFlipped; 

foi = 20; % ssvef freq 
toi = abs(p.tstart) + p.eventTimes(1):abs(p.tstart) + p.eventTimes(4); % precue to response cue 

subjectIdx = 1:numel(subjectNames); 
% subjectIdx(5) = []; % if doing 9 subjects 
nSubs = numel(subjectIdx); 

%% 
data1 = squeeze(data1(foi,toi,subjectIdx)); 
data2 = squeeze(data2(foi,toi,subjectIdx)); 

