% meg_runGroup 

%% setup
user = 'mcq'; % 'mcq','karen'
expt = 'TA2';

%% get session names
switch expt
    case 'TANoise'
        sessionNames = {'R0817_20171212','R0817_20171213',...
            'R1187_20180105','R1187_20180108',...
            'R0983_20180111','R0983_20180112',...
            'R0898_20180112','R0898_20180116',...
            'R1021_20180208','R1021_20180212',...
            'R1103_20180213','R1103_20180215',...
            'R0959_20180219','R0959_20180306',...
        	'R1373_20190723','R1373_20190725',...
            'R1452_20190717','R1452_20190718',...
            'R1507_20190702','R1507_20190705'}; % N=10 x 2 sessions TANoise

    case 'TA2'
        sessionNames = {'R0817_20181120', 'R0817_20190625',...
            'R0898_20190723', 'R0898_20190724',...
            'R0959_20181128', 'R0959_20190703',...
            'R0983_20190722', 'R0983_20190723',...
            'R1103_20181121', 'R1103_20190710',...
            'R1187_20181119', 'R1187_20190703',...
            'R1373_20181128', 'R1373_20190708',...
            'R1452_20181119', 'R1452_20190711',...
            'R1507_20190621', 'R1507_20190627',...
            'R1547_20190729', 'R1547_20190730'}; % N=10 x 2 sessions TA2
end

%%% RD suggestion: Instead of having many blocks of sessionNames, or 
%%% deleting sessions later, select from the above lists using indices
sessions = 1:16;
sessionNames = sessionNames(sessions);

%% TA2 March 5 block mat prep 
sessionNames = {'R0898_20190723','R0898_20190724',...
    'R0959_20181128','R0959_20190703',...
    'R0983_20190722','R0983_20190723',...
    'R1103_20181121','R1103_20190710',...
    'R1187_20181119','R1187_20190703',...
    'R1373_20181128','R1373_20190708',...
    'R1452_20181119','R1452_20190711',...
    'R1547_20190729','R1547_20190730'} 
%     'R1507_20190621','R1507_20190627',...
%     'R1535_20190708','R1535_20190711',...
%     ,...

sessionNames = {'R1507_20190621','R1507_20190627','R1535_20190708','R1535_20190711'};

%% TA2 minus 1
sessionNames = {'R0817_20181120', 'R0817_20190625',...
    'R0898_20190723','R0898_20190724',...
    'R0959_20181128','R0959_20190703',...
    'R0983_20190722','R0983_20190723',...
    'R1103_20181121','R1103_20190710',...
    'R1187_20181119','R1187_20190703',...
    'R1373_20181128','R1373_20190708',...
    'R1452_20181119','R1452_20190711',...
    'R1507_20190621','R1507_20190627',...
    'R1535_20190708','R1535_20190711',...
    'R1547_20190729','R1547_20190730',...
    };

%% TA2 batch 1
% TA2
sessionNames = {'R0817_20181120', 'R0817_20190625',...
    'R0890_20181121',...
    'R0898_20190723','R0898_20190724',...
    'R0959_20181128','R0959_20190703',...
    'R0983_20190722','R0983_20190723',...
    'R1103_20181121','R1103_20190710',...
    'R1187_20181119','R1187_20190703',...
    'R1373_20181128','R1373_20190708',...
    'R1452_20181119','R1452_20190711',...
    'R1507_20190621','R1507_20190627',...
    'R1535_20190708','R1535_20190711',...
    'R1547_20190729','R1547_20190730',...
    };


%% TA2 batch 1 
sessionNames = {'R0817_20181120', 'R0817_20190625',...
    'R0890_20181121',...
    'R0898_20190723','R0898_20190724',...
    'R0959_20181128','R0959_20190703',...
    'R0983_20190722','R0983_20190723',...
    'R1103_20181121'}; 

%% TA2 tiny batch check 
sessionNames = {'R0817_20181120','R0817_20190625'}; 
   
%% TA2 batch 2
sessionNames = {'R1103_20190710',...
    'R1187_20181119','R1187_20190703',...
    'R1373_20181128','R1373_20190708',...
    'R1452_20181119','R1452_20190711',...
    'R1507_20190621','R1507_20190627',...
    'R1535_20190708','R1535_20190711',...
    'R1547_20190729','R1547_20190730',...
    };

%% TANoise 
sessionNames = {'R0817_20171212', 'R0817_20171213',...
    'R0898_20180112','R0898_20180116',...
    'R0959_20180219','R0959_20180306',...
    'R0983_20180111','R0983_20180112',...
    'R1021_20180208','R1021_20180212',...
    'R1103_20180213','R1103_20180215',...
    'R1187_20180105','R1187_20180108',...
    'R1373_20190723','R1373_20190725',...
    'R1452_20190717','R1452_20190718',...
    'R1507_20190702','R1507_20190705'
    };

%% run analysis 
makeGroup = 0; 
 
groupA = [];
for i=1:numel(sessionNames)
    sessionDir = sessionNames{i}; 
    disp(sessionDir)

    [A, D, selectedChannels] = meg_runAnalysis(expt, sessionDir, user); 
    groupA{i} = A;
    pause(.1)
%     close all
    
    if makeGroup
        groupD(i).sessionDir = sessionDir;
        groupD(i).data = D; % data by cue cond
        groupD(i).selectedChannels = selectedChannels;
    end
end


%% alpha finder
for i=1:numel(sessionNames)
    sessionDir = sessionNames{i}; 
    disp(sessionDir)
    meg_selectChannels(sessionDir)
    close all
end

%% delete sessions

sessionsDeleteIdx = [3,22,23]; 
sessionsDelete = sessionNames(sessionsDeleteIdx); 
groupD(sessionsDeleteIdx) = []; % clear rejected sessions 
sessionNames(sessionsDeleteIdx) = []; % clear rejected session names 

%% make group data structure 2

% how many channels to inspect? (max 7, if >7 rerun analyses) 
nC = 5; 

% cue conds 
groupD2.cueT1 = []; 
groupD2.cueT2 = []; 
groupD2.neutral  = []; 
groupD2.cueT1All = []; 
groupD2.cueT2All = []; 
groupD2.neutralAll = []; 

fields = fieldnames(groupD2); 

% create structure groupD2 of group cond data
for i=1:numel(sessionNames)
    for f=1:numel(fields)
        vals = groupD(i).data.(fields{f});
        selectedChannelsVals = vals(:,groupD(i).selectedChannels(1:nC),:);  % time x channels x trials
        meanVals = nanmean(selectedChannelsVals,3); % mean across time
        groupD2.(fields{f}) = cat(3,groupD2.(fields{f}),meanVals);
    end
end

% groupD2 all time x channel x trial x subject
for i=1:numel(sessionNames)
    for f=1:3 
        vals = groupD(i).data.(fields{f});
        selectedChannelsVals = vals(:,groupD(i).selectedChannels(1:nC),:);  % time x 7 channels x trials
        groupD2.([fields{f} 'All']) = cat(4,groupD2.([fields{f} 'All']),selectedChannelsVals); 
    end
end

% average group
groupD2.meanCueT1    = nanmean(nanmean(groupD2.cueT1,3),2);
groupD2.meanCueT2    = nanmean(nanmean(groupD2.cueT2,3),2);
groupD2.meanNeutral  = nanmean(nanmean(groupD2.neutral,3),2);

% group difference calcs 
groupD2.cueT1_neutral = groupD2.meanCueT1 - groupD2.meanNeutral; 
groupD2.cueT2_neutral = groupD2.meanCueT2 - groupD2.meanNeutral; 
groupD2.cueT1_cueT2 = groupD2.meanCueT1 - groupD2.meanCueT2; 
groupD2.cueT2_cueT1 = groupD2.meanCueT2 - groupD2.meanCueT1; 

save('groupD2.mat','groupD2')

%% make subject structure (avg across sessions) 

nSubjects = numel(sessionNames)/2; 

subjectD2.cueT1 = []; 
subjectD2.cueT2 = []; 
subjectD2.neutral  = []; 

subjectD2.cueT1All = []; 
subjectD2.cueT2All = []; 
subjectD2.neutralAll = []; 

subjectD2.cueT1Trial = []; 
subjectD2.cueT2Trial = []; 
subjectD2.neutralTrial = []; 

fields = fieldnames(subjectD2); 

for i = 1:2:20
    for f=13:15 % 3 cue conds, all trials (fields 1:3 for cue conds avg) 
        subjectVals = nanmean(groupD2.(fields{f})(:,:,:,i:i+1),4);
        subjectD2.(fields{f}) = cat(4,subjectD2.(fields{f}),subjectVals);
    end
end


%% 
% average across channels 
subjectD2.meanCueT1    = squeeze(nanmean(subjectD2.cueT1,2));
subjectD2.meanCueT2    = squeeze(nanmean(subjectD2.cueT2,2));
subjectD2.meanNeutral  = squeeze(nanmean(subjectD2.neutral,2));

%% setup group figures 

p = meg_params('TA2_Analysis'); 
t = p.tstart:p.tstop;
xlims = [p.tstart,p.tstop]; 
colorDiff = [0.8 0.6 1.0]; 

%% group figures 
%% ERF avg 

fields = fieldnames(groupD2); 
figure % ERF by group 
set(gcf,'Position',[100 100 1200 500])
hold on
for f = 4:6 % 3 cue conds mean
    l = plot(t,groupD2.(fields{f}),'LineWidth',2);
    l.Color = [p.cueColors(f-3,:) p.colorAlpha]; 
end
yline(0,'-k');
xlim(xlims)
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT1','cueT2','neutral')
title(sprintf('TA2 group top %d visual channels (by avg peak prom) ERF',nC))
vline(p.eventTimes,'k',p.eventNames)

%% ERF group by session 

figure
figpos = 1; 
for i = 1:numel(sessionNames)
    subplot(numel(sessionNames)/2,2,figpos)
    hold on
    for f = 1:3 % 3 cue conds 
        l = plot(t,nanmean(groupD2.(fields{f})(:,:,i),2),'LineWidth',2);
        l.Color = [p.cueColors(f,:) p.colorAlpha];
    end
    yline(0,'-k');
    xlim(xlims)
    vline(p.eventTimes,'k')
    title(sprintf('%s',und2space(sessionNames{i})))
    figpos = figpos + 1;
end
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT1','cueT2','neutral')
rd_supertitle2(sprintf('TA2 group top %d visual channels (by avg peak prom) ERF',nC))

%% ERF by subject, avg across 2 sessions 

figure
for i = 1:nSubjects
    subplot(nSubjects,1,i)
    hold on
    for f = 4:6
        l = plot(t,subjectD2.(fields{f})(:,i)); 
        l.Color = [p.cueColors(f-3,:) p.colorAlpha];
        
         % peak stuff 
         [pks,locs,w,prom] = findpeaks(subjectD2.(fields{f})(:,i));
         [sortProm, promIdx] = sort(prom,'descend');
         plot(promIdx(1:2),sortProm(1:2),'.b', 'MarkerSize', 10)
    end
    yline(0,'-k');
    xlim(xlims)
    vline(p.eventTimes,'k')
    title(sprintf('%s',und2space(sessionNames{i})))
end
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT1','cueT2','neutral')
rd_supertitle2(sprintf('TA2 avg ERF by subject top %d visual channels',nC))

%% difference figs

figure
subplot(4,1,1)
hold on
l = plot(t,groupD2.cueT1_neutral,'LineWidth',2); 
l.Color = [0.6 0.8 1]; 
yline(0,'-k');
xlim(xlims)
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT1-neutral')
vline(p.eventTimes,'k',p.eventNames)

subplot(4,1,2)
hold on
l = plot(t,groupD2.cueT2_neutral,'LineWidth',2); 
l.Color = [1 0.6 0.6]; 
yline(0,'-k');
xlim(xlims)
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT2-neutral')
vline(p.eventTimes,'k',p.eventNames)

subplot(4,1,3)
hold on
l = plot(t,groupD2.cueT1_cueT2,'LineWidth',2);
l.Color = [0.4940, 0.1840, 0.5560 0.7] 
yline(0,'-k');
xlim(xlims)
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT1-cueT2')
vline(p.eventTimes,'k',p.eventNames)

subplot(4,1,4)
hold on
l = plot(t,groupD2.cueT2_cueT1,'LineWidth',2); 
l.Color = [0.4940, 0.1840, 0.5560 0.7] 
yline(0,'-k');
xlim(xlims)
xlabel('time (ms)')
ylabel('amp (T)')
legend('cueT2-cueT1')
vline(p.eventTimes,'k',p.eventNames)

rd_supertitle2('TA2 group top 7 visual channels (by avg peak prom) ERF')

%% peaks by subject by cue 

for f = 14:16 % 3 cue cond 
    [pks,locs,w,prom] = findpeaks(subjectD2.(field{f}));
    [sortProm, promIdx] = sort(prom,'descend'); 
end

%% max T1 window by subject by cue 

% windowT1 = [p.eventTimes(2):p.eventTimes(2)+300]; 
idxT1 = find(t==p.eventTimes(2));
window = idxT1:idxT1+300; 

subjectD2.maxCueT1 = [];
subjectD2.maxCueT2 = [];
subjectD2.maxNeutral = []; 

fields = fieldnames(subjectD2); 

%% max by subject
figure
for i = 1:nSubjects
    subplot(nSubjects,1,i)
    hold on
    for f = 4:6
        l(f-3) = plot(t,subjectD2.(fields{f})(:,i)); 
        l(f-3).Color = [p.cueColors(f-3,:) p.colorAlpha];
        hold on
        % max
        [M,I] = max(subjectD2.(fields{f})(window,i));
        plot(t(window(I)),M,'.b','MarkerSize', 10)
        
        maxVal = [M; t(window(I))]; 
        subjectD2.(fields{f+3}) = cat(2,subjectD2.(fields{f+3}),maxVal);
    end
    plot(t(window),zeros(size(t(window))),'k','LineWidth',2) % plot window time
    yline(0,'-k');
    xlim(xlims)
    vline(p.eventTimes,'k')
    title(sprintf('%s',und2space(sessionNames{i})))
end
xlabel('time (ms)')
ylabel('amp (T)')
legend(l,{'cueT1','cueT2','neutral'})
rd_supertitle2(sprintf('TA2 avg ERF by subject top %d visual channels',nC))

%% check T1 timing and max cueCond

figure
for i = 1:nSubjects
    subplot(1,nSubjects,i) 
    hold on 
    maxN = subjectD2.maxNeutral(1,i); 
    scatter(subjectD2.maxCueT1(2,i),subjectD2.maxCueT1(1,i)/maxN,30,p.cueColors(1,:),'filled')
    scatter(subjectD2.maxCueT2(2,i),subjectD2.maxCueT2(1,i)/maxN,30,p.cueColors(2,:),'filled')
    scatter(subjectD2.maxNeutral(2,i),maxN/maxN,30,p.cueColors(3,:),'filled')
    xlim(xlims)
    ylim([0.8 1.25])
end
rd_supertitle2('peakT1')


%% move files
% sourceDir = '/e/1.3/p1/denison/Downloads/MEG';
% exptDir = meg_pathToTAMEG(expt, user);
% 
% success = [];
% for i=1:numel(sessionNames)
%     if success(i)==0
%         sessionDir = sessionNames{i};
%         disp(sessionDir)
%         
%         dataDir = sprintf('%s/%s', exptDir, sessionDir);
%         matDir = sprintf('%s/mat', dataDir);
%         
%         fileBase = meg_sessionDirToFileBase(sessionDir, expt);
%         dataFile = sprintf('%s_bietfp_data.mat', fileBase);
%         
%         sourceFile = sprintf('%s/%s', sourceDir, dataFile);
%         destFile = sprintf('%s/%s', matDir, dataFile);
%         
%         tic
%         [success(i), message{i}] = movefile(sourceFile, destFile);
%         toc
%     end
% end
