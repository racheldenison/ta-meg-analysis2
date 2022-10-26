sessionNames = meg_sessions('TANoise'); 

for i = 11:numel(sessionNames)
    disp(sessionNames{i})
    meg_runPreproc(sessionNames{i}) 
end