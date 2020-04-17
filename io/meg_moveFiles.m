% meg_moveFiles.m

%% move files
sourceDir = '/e/1.3/p1/denison/Downloads/MEG';
exptDir = meg_pathToTAMEG(expt, user);

success = zeros(1,numel(sessionNames));
for i=1:numel(sessionNames)
    if success(i)==0
        sessionDir = sessionNames{i};
        disp(sessionDir)
        
        dataDir = sprintf('%s/%s', exptDir, sessionDir);
        matDir = sprintf('%s/mat', dataDir);
        
        fileBase = meg_sessionDirToFileBase(sessionDir, expt);
        dataFile = sprintf('%s_bietfp_data.mat', fileBase);
        
        sourceFile = sprintf('%s/%s', sourceDir, dataFile);
        destFile = sprintf('%s/%s', matDir, dataFile);
        
        tic
        [success(i), message{i}] = movefile(sourceFile, destFile);
        toc
    end
end

%% save individual channel selection files
exptDir = meg_pathToTAMEG('TA2');
load(sprintf('%s/Group/mat/groupC.mat',exptDir))
load vis/data_hdr.mat

for i = 1:numel(groupC)
    C = groupC(i);
    C.channelsRanked = C.sortChByProm; % just make a clearer fieldname
    sessionDir = C.sessionDir;
    disp(sessionDir)
    
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    fileName = sprintf('%s/channels_peakprom.mat',matDir);
%     save(fileName,'C')

    [Y, idx] = sort(C.channelsRanked);
    ssm_plotOnMesh(C.sortProm(idx)', '', [], data_hdr, '2d');
    pause(1)
end