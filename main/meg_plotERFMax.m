function meg_plotERFMax 

% finds max amp and max idx per target per condition per subject
% averages and plots 
% April 2020 

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

targetNames = p.eventNames(2:3); 
fields = fieldnames(D.subject); 
windowSize = 300;  % window size to look for target evoked max 

[sessionNames,subjectNames] = meg_sessions(expt); 


%% setup A structure 
%  stores max amp and time
%  by A.targetName.cueName: trials x subjects x max
A = []; M = []; I = []; 
for iT = 1:nTargets
    for iF = 1:numel(fields) 
        A.(targetNames{iT}).(cueConds{iF}) = [];
    end
end

%% find max and max idx per target per condition per subject 

for iF = 1:numel(fields) % cue condition
    for iT = 1:numel(targetNames) % target
        for iS = 1:numel(subjectNames) % subject
            target = targetNames(iT);
            idxTargetName = find(strcmp(p.eventNames,target));
            idxT = find(t==p.eventTimes(idxTargetName));
            window = idxT:idxT+windowSize;
            
            % by subject
            [M(iS),I(iS)] = max(D.subject.(fields{iF})(window,iS));
            A.(targetNames{iT}).(fields{iF}).max = M;
            A.(targetNames{iT}).(fields{iF}).maxIdx = t(window(I));
            
        end
        % group avg
        [groupM,groupI] = max(D.group.(fields{iF})(window));
        A.(targetNames{iT}).(fields{iF}).groupMax = groupM;
        A.(targetNames{iT}).(fields{iF}).groupMaxIdx = t(window(groupI));
    end
end

% compute max and max idx averages and std across subjects 
for iT = 1:numel(targetNames)
    for iF = 1:numel(fields)
        % max 
        A.(targetNames{iT}).(fields{iF}).maxMean = nanmean(A.(targetNames{iT}).(fields{iF}).max,2);
        A.(targetNames{iT}).(fields{iF}).maxStd = nanstd(A.(targetNames{iT}).(fields{iF}).max,[],2);
        % max idx 
        A.(targetNames{iT}).(fields{iF}).maxIdxMean = nanmean(A.(targetNames{iT}).(fields{iF}).maxIdx,2);
        A.(targetNames{iT}).(fields{iF}).maxIdxStd = nanstd(A.(targetNames{iT}).(fields{iF}).maxIdx,[],2);
    end
end

%% scatterplot max w v and h error bars 

for iT = 1:numel(targetNames)
    target = targetNames(iT);
    idxTargetName = find(strcmp(p.eventNames,target));
    idxT = find(t==p.eventTimes(idxTargetName));
    window = idxT:idxT+windowSize;
    
    figure
    set(gcf,'Position',[100 100 500 500])
    hold on
    for iF = 1:numel(fields)
        % vertical error bar
        errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
            A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
            A.(targetNames{iT}).(fields{iF}).groupMax/numel(subjectNames),... % vertical error bar: max ste
            'Color',p.cueColors(iF,:),'LineWidth',2);
        % horizontal error bar
        errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
            A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
            A.(targetNames{iT}).(fields{iF}).groupMaxIdx/numel(subjectNames),... % vertical error bar: max ste
            'horizontal','Color',p.cueColors(iF,:),'LineWidth',2);
        
        p1 = plot(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,...
            A.(targetNames{iT}).(fields{iF}).groupMax,...
            'o','MarkerSize',15);
        p1.Color = p.cueColors(iF,:);
        p1.MarkerFaceColor = p.cueColors(iF,:);
    end
    ylim([0.5e-13,1.8e-13])
    xlim([t(window(1))-50,t(window(end))])
    xlabel('time (ms')
    ylabel('amp (T)') 
    vline(p.eventTimes(iT+1),'k:','target') 
    title(sprintf('%s max by cue condition',targetNames{iT}))
end

%% plot ERF w max scatter 

p.colorTarget = [0.5 0.5 0.5; 0 0 0];

figure
set(gcf,'Position',[100 100 1200 350])
hold on
for iF = numel(fields):-1:1
    % plot ERF
    p1 = plot(t,D.group.(fields{iF}),'LineWidth',2);
    p1.Color = p.cueColors(iF,:);
    p1.Color(4) = p.colorAlpha;
    for iT = 1:numel(targetNames)
        target = targetNames(iT);
        idxTargetName = find(strcmp(p.eventNames,target));
        idxT = find(t==p.eventTimes(idxTargetName));
        window = idxT:idxT+windowSize;
        % vertical error bar
        errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
            A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
            A.(targetNames{iT}).(fields{iF}).groupMax/numel(subjectNames),... % vertical error bar: max ste
            'Color',p.cueColors(iF,:),'LineWidth',0.5);
        % horizontal error bar
        errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
            A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
            A.(targetNames{iT}).(fields{iF}).groupMaxIdx/numel(subjectNames),... % vertical error bar: max ste
            'horizontal','Color',p.cueColors(iF,:),'LineWidth',0.5);
        plot(t(window),zeros(size(window)),'LineWidth',4,'Color',p.colorTarget(iT,:))
    end
end
hline(0,'k--')
xlim([t(1),t(end)])
xlabel('time (ms')
ylabel('amp (T)')
vline(p.eventTimes,'k--',p.eventNames)
title(sprintf('group max by target by cue condition'))

