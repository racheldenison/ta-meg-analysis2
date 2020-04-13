% meg_runAnalysisMultiSession

%% Setup
expt = 'TA2';
% sessionDirs = {'R0817_20181120','R0817_20190625'};
sessionDirs = {'R1187_20181119','R1187_20190703'};
user = 'mcq';
paramType = 'Analysis';

% params
p = meg_params(sprintf('%s_%s', expt, paramType));

nSessions = numel(sessionDirs);

saveAnalysis = 0;
saveFigs = 0;

%% Get individual session data structures
clear D I B
for iS = 1:nSessions
    sessionDir = sessionDirs{iS};
    disp(sessionDir)
    
    [~, Ds(iS), ~, Is(iS), Bs(iS)] = meg_runAnalysis(expt, sessionDir, user);
end

%% Make multi session data structures
dFields = fieldnames(Ds(1));
bFields = fieldnames(Bs(1));

clear D I B
% initialize super structures
for iF = 1:numel(dFields)
    fname = dFields{iF};
    D.(fname) = [];
    I.(fname) = [];
end

B.responseData_labels = Bs(1).responseData_labels;
for iF = 1:numel(bFields)
    fname = bFields{iF};
    if isnumeric(Bs(1).(fname))
        B.(fname) = [];
    end
end

% aggregate sessions
for iS = 1:nSessions
    for iF = 1:numel(dFields)
        fname = dFields{iF};
        D.(fname) = cat(3, D.(fname), Ds(iS).(fname));
        I.(fname) = cat(1, I.(fname), Is(iS).(fname));
    end
    
    for iF = 1:numel(bFields)
        fname = bFields{iF};
        if isnumeric(Bs(1).(fname))
            B.(fname) = cat(1, B.(fname), Bs(iS).(fname));
        end
    end
end

%% Decoding orientation
analysisName = 'classAcc';
targetNames = {'T1','T2'};
classNames = {'V','H'};
twin = [-50 550];
condNames = fields(D);
selectedChannels = [];

for iT = 1:2
    targetTime = p.eventTimes(strcmp(p.eventNames,targetNames{iT}));
    p.targetWindow = targetTime + twin;
    classLabels = B.t1t2Axes(:,iT);
    [A(iT), fH{iT}, figNames{iT}] = meg_plotDecode(D, I, p, classLabels, classNames, selectedChannels);
    
    if saveFigs
        for iF = 1:numel(fH{iT})
            figNamesT{iF} = sprintf('%s_%s', figNames{iT}{iF}, targetNames{iT});
        end
        rd_saveAllFigs(fH{iT}, figNamesT, [], figDir)
    end
end

if saveAnalysis
    decodeAnalStr = A(1).(condNames{1}).decodingOps.analStr;
    %             save(sprintf('%s/%s_%s_%sSlice_%s.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'A', 'targetNames')
    %             save(sprintf('%s/%s_%s_%sSlice_%s_varyNChannels10-157_sortedBy20Hz.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr), 'allA', 'targetNames', 'nTopCh')
    save(sprintf('%s/%s_%s_%sSlice_%s_nCh%d.mat', matDir, analysisName, analStr, sliceType, decodeAnalStr, nTopCh(i)), 'A', 'targetNames')
end
