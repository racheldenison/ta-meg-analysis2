function meg_plotGroupMax 

load('/Users/kantian/Dropbox/github/ta-meg-analysis2/groupDataMat/TA2/D.mat')

% load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/mat/subjectD2.mat')
% 
% D.channelSortMethod = 'peakProm'; 
% D.nChannels = 5; 
% 
% D.subject.cueT1 = subjectD2.meanCueT1; 
% D.subject.cueT2 = subjectD2.meanCueT2; 
% D.subject.cueN  = subjectD2.meanNeutral; 
% 
% D.group.cueT1   = subjectD2.groupCueT1; 
% D.group.cueT2   = subjectD2.groupCueT2;  
% D.group.cueN    = subjectD2.groupNeutral;  

%% setup 

expt = 'TA2'; 
p = meg_params('TA2_Analysis'); 

tIdx = 1:p.trialTime; % epoch time, 1 indexed
t = p.tstart:p.tstop; % epoch time, relative to precue

nTargets = 2; 
targetNames = p.eventNames(2:3); 
cueConds = fieldnames(D.subject); 
windowSize = 300;  % window size to look for target evoked max 

[sessionNames,subjectNames] = meg_sessions(expt); 
nSessions = numel(sessionNames); 
nSubjects = numel(subjectNames); 

%% setup A structure 
%  stores max amp and time by A.targetName.cueName
A = []; 
for iT = 1:nTargets
    for c = 1:numel(cueConds) 
        A.(targetNames{iT}).(cueConds{c}) = [];
    end
end

%% max 
%% window by subject by cue 

for iT = 1:nTargets 
    target = targetNames(iT); 
    idxTargetName = find(strcmp(p.eventNames,target)); 
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowSize; 
end


%% max T1 

fields = {'meanCueT1','meanCueT2'}; 

for s = 1:10
    maxMat = [];
        [M,I] = max(A.(fields{nField})(window,s));
        % plot(t(window(I)),M,'.b','MarkerSize', 10)
        maxVal = [M, iT(window(I))];
        maxMat = cat(1,maxMat,maxVal);
    subMax(s,:,:) = maxMat; 
end
% avgMax = squeeze(nanmean(subMax,2)); 
% avgMax = squeeze(avgMax(:,1)); 

I = subMax(:,:,2); 
avgI = squeeze(nanmean(I,1)); 
stdI = nanstd(I,0,1); 
steI = stdI/sqrt(numel(subjectNames)); 

M = subMax(:,:,1); 
avgM = squeeze(nanmean(M,1)); 
stdM = nanstd(M,0,1); 
steM = stdM/sqrt(numel(subjectNames)); 

% store into maxT1 structure (cueT1)
maxT1Subject.cueT1.I = I;
maxT1Subject.cueT1.avgI = avgI; 
maxT1Subject.cueT1.stdI = stdI;
maxT1Subject.cueT1.steI = steI;

maxT1Subject.cueT1.M = M;
maxT1Subject.cueT1.avgM = avgM; 
maxT1Subject.cueT1.stdM = stdM;
maxT1Subject.cueT1.steM = steM;

% store into maxT1 structure (cueT2)
maxT1Subject.cueT2.I = I;
maxT1Subject.cueT2.avgI = avgI; 
maxT1Subject.cueT2.stdI = stdI;
maxT1Subject.cueT2.steI = steI;

maxT1Subject.cueT2.M = M;
maxT1Subject.cueT2.avgM = avgM; 
maxT1Subject.cueT2.stdM = stdM;
maxT1Subject.cueT2.steM = steM;

% % store into maxT1 structure (cueNeutral)
% maxT1Subject.cueNeutral.I = I;
% maxT1Subject.cueNeutral.avgI = avgI; 
% maxT1Subject.cueNeutral.stdI = stdI;
% maxT1Subject.cueNeutral.steI = steI;
% 
% maxT1Subject.cueNeutral.M = M;
% maxT1Subject.cueNeutral.avgM = avgM; 
% maxT1Subject.cueNeutral.stdM = stdM;
% maxT1Subject.cueNeutral.steM = steM;

%% plot setup

window = p.eventTimes(2):p.eventTimes(2)+300; 

%% plot

fields = {'cueT1','cueT2'}; 
nSubjects = numel(subjectNames); 

figure
set(gcf,'Position',[100 100 2000 500]) 
hold on 
for s = 1:nSubjects
    subplot(1,nSubjects,s)
    xlim([window(1),window(end)])
    for f = 1:numel(fields)
        condName = fields{f}; 
        condVals = maxT1Subject.(fields{f}); 
        
        hold on
        errV = errorbar(condVals.avgI(s),condVals.avgM(s),condVals.steM(s),'Color',p.cueColors(f,:),'LineWidth',2) % vertical error bar
        errH = errorbar(condVals.avgI(s),condVals.avgM(s),condVals.steI(s),'horizontal','Color',p.cueColors(f,:),'LineWidth',2) % horizontal error bar
        
        errV.Color(4) = p.colorAlpha;
        errH.Color(4) = p.colorAlpha;
         
        p1 = plot(condVals.avgI(s),condVals.avgM(s),'MarkerFaceColor',p.cueColors(f,:))
        p1.Color(4) = p.colorAlpha; 
    end
    title(sprintf('%s',subjectNames{s}))
    ylabel('max T1 amp (T)')
    xlabel('epoch time (ms)')
    % legend(p1,fields) fix this 
    rd_supertitle2('max amp in T1 window')
end

%% plot 

A.subMax = subMax; 

figure
hold on
plot(avgI,avgM,'o') 
errorbar(avgI,avgM,steM,'o')
errorbar(avgI,avgM,steI,'horizontal','o')



