function subjectData = meg_sessions2subjects(sessionData)
% function subjectData = meg_sessions2subjects(sessionData)
% takes 1 x 20 (session) data and averages to 1 x 10 (subject) level 

dims = size(sessionData); 
nSessions = dims(2); 
if numel(dims) >2 
    error('Check size of input. Should be val x 20 sessions.')
end
nSessionsPerSubject = 2; % number of sessions per subject 

sCount = 1;
for i = 1:nSessionsPerSubject:nSessions
    clear val
        sVal = sessionData(:,i:i+nSessionsPerSubject-1);
        sVal = mean(sVal,2,'omitnan');
        subjectData(:,sCount) = sVal;
    sCount=sCount+1;
end

