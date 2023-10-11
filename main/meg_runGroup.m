% meg_runGroup 
% group meg_runAnalysis 

% add ft 
% addpath '/Users/kantian/Google Drive/Software/fieldtrip-20190416'
% ft_defaults

%% setup
user = 'karen'; % 'mcq','karen','rachel'
expt = 'TA2'; % ' TANoise' 'TA2'
sessionIdxs = 1:2; 
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

%% get session names
% allSessions = meg_sessions(expt);
% sessionNames = allSessions(sessionIdx);

%% run analysis 
clear groupA
clear groupD
clear groupI
clear groupB

%%
for i = sessionIdxs % :20 %20 15:numel(sessionNames) % 8:numel(sessionNames) % 1:numel(sessionNames), skipped 8, cant read sqd?
    sessionDir = sessionNames{i}; 
    disp(sessionDir)
    sessionIdx = i; 
    
    [A, D, selectedChannels, I, B] = meg_runAnalysis(expt,sessionDir,user); 
    % [A] = loadGroupAnalysis(expt,sessionDir,user); 
    %     A2 = [];
    %     A2.tfAmps = A.all.tfAmps;
    %     A2.tfPows = A.all.tfPows;
    % groupITPC_40Hz_avgTrial{i} = A;
    % group_TF_wholeTrial{i} = A;
    % group_TF_preTrial{i} = A;
    % group_TF_preTarget{i} = A;
    %     groupD{i} = D;
    % groupB{i} = B;
    % groupITPCspectrogram{i} = A.all.ITPCMean; 
    % groupTFspectrogram_avgTrial(i) = A; 
%     groupTFspectrogram_singleTrial(i).normAmps = A.all.normAmps; 
%     groupTFspectrogram_singleTrial(i).normPows = A.all.normPows; 

%     groupA(i).cueT1.normAmps = A.cueT1.normAmps; 
%     groupA(i).cueT1.normPows = A.cueT1.normPows; 
%     groupA(i).cueT1.amps = A.cueT1.meanTfAmps; 
%     groupA(i).cueT1.pows = A.cueT1.meanTfPows; 
%     
%     groupA(i).cueT2.normAmps = A.cueT2.normAmps; 
%     groupA(i).cueT2.normPows = A.cueT2.normPows; 
%     groupA(i).cueT2.amps = A.cueT2.meanTfAmps; 
%     groupA(i).cueT2.pows = A.cueT2.meanTfPows; 

%     groupA(i).cueT1.ITPCMean = A.cueT1.ITPCMean; 
%     groupA(i).cueT2.ITPCMean = A.cueT2.ITPCMean;  
    
    % groupB(i) = B; 
    groupA(i) = A; 
    % groupI(i) = I; 
    % groupD(i) = D; 
    close all 
end
% groupA.type = 'ITPCspectrogram_singleTrial'; 

%% combine 
valsAll = []; 
for i = 1:numel(sessionNames)
    % vals = groupA(i).cueT2.ITPCMean; 
    vals = groupTFspectrogram_avgTrial(i).cueT2.pows;
    valsAll = cat(3,valsAll,vals); 
end
vals_cueT2 = valsAll; 

%% combine behavior 
for i = 1:numel(sessionNames)
    discrimHMFC(:,:,i) = groupB(i).discrimHMFC; 
    t1t2Axes(:,:,i) = groupB(i).t1t2Axes; 
end

%% behavior d' 
for i = 1:numel(sessionNames)
    subjectDprime(i) = norminv(discrimHMFC(:,1,i)); 
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

%% ITPC breakdown into subjects 
slice = 'all'; 
sessionVals = []; 
for i = 1:numel(sessionNames)
    vals = groupA(i).(slice).ITPCMean; 
    sessionVals(:,:,i) = vals; 
end
A.all.session = []; 
A.all.session = sessionVals; 

%% average session to subjects 
variableOld = 'normSession'; 
variableNew = 'normSubject'; 
conds = {'cueT1','cueT2','all'}; 

count = 1; 
for i=1:2:numel(sessionNames)
    for iC = 1:numel(conds)
        subjectVals = [];
        vals = A.(conds{iC}).(variableOld)(:,:,i:i+1);
        vals = mean(vals,3,'omitnan');
        subjectVals = cat(3,subjectVals,vals);
        A.(conds{iC}).(variableNew)(:,:,count) = subjectVals;
    end
    count = count+1; 
end

%% normalize 
% find and store baseline mean per frequency 
for i = 1:numel(sessionNames)
    baselineVals = A.all.session(:,p.baselineNorm+abs(p.tstart),i); 
    baselineVals = mean(baselineVals,2,'omitnan'); % average baseline across time 
    A.all.baselineITPC(:,i) = baselineVals; 
end

% subtract baseline for all trials 
slice = 'all'; 
conds = {'cueT1','cueT2'}; 
variable = 'normSession';
A.cueT1.normSession = []; 
A.cueT2.normSession = []; 
for i = 1:numel(sessionNames)
    A.(slice).(variable)(:,:,i) = A.(slice).session(:,:,i)-A.(slice).baselineITPC(:,i); 
    for iC = 1:2
        A.(conds{iC}).(variable)(:,:,i) = A.(conds{iC}).session(:,:,i)-A.(slice).baselineITPC(:,i); 
    end
end
% A.(slice).normSession
% A = rmfield(A.cueT1,'baselineVal'); 

%% flip subjects 
for i = 1:10
    for iC = 1:numel(conds)
        if ITPCsubject(i) == -1 % flip if downer
            A.(conds{iC}).normSubjectFlipped(:,:,i) = A.(conds{iC}).normSubject(:,:,i)*-1; 
        else % otherwise just copy 
            A.(conds{iC}).normSubjectFlipped(:,:,i) = A.(conds{iC}).normSubject(:,:,i); 
        end
    end
end

%% save
% save('groupA_ITPCspectrogram_byAtt.mat','A','-v7.3')

