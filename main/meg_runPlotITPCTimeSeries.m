% meg_runPlotITPCTimeSeries 



for i = 20:-1:1
    sessionDir = sessionNames{i}; 
    [fH,figNames] = meg_plotITPCTimeSeries(sessionDir,i); 
    close all
end