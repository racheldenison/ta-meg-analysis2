function subjectData = meg_sessions2subjects(sessionData)
% function subjectData = meg_sessions2subjects(sessionData)
% takes 1 x 20 (session) data and averages to 1 x 10 (subject) level 

nSessions = size(sessionData,2); 
nSessionsPerSubject = 2; % number of sessions per subject 

sCount = 1; 
for i = 1:nSessionsPerSubject:nSessions
    clear sVal
    sVal = sessionData(i:i+nSessionsPerSubject-1);
    sVal = mean(sVal);
    subjectData(sCount) = sVal; 
    sCount=sCount+1; 
end
