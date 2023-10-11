function [D,A] = meg_peakAnalysis(D4,nChOI,groupC,expt,user,slice)

% function [D,A] = meg_peakAnalysis(D4,nChOI,groupC,expt,user,slice)
% Finds peak amplitudes and timings of targets evoked response by condition
% averages channels (nChOI) > per condition per subject
% averages and plots 
% April 2020 

% D4 = data structure, fields = slice conditions 
%   D4.cue (time? x channels x session) 

% mega group D 
% load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/mat/groupD.mat') 
% load('/Users/kantian/Dropbox/github/ta-meg-analysis2/groupDataMat/TA2/D.mat')
% load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/mat/subjectD2.mat')
%% checks 
if ~exist('user','var')
    user = [];
end
if ~exist('expt','var')
    error('expt type (TA2 or Noise) must be defined')
end
if ~exist('slice','var')
    slice = 'cue'; 
    disp('slice not specified, default by cue')
end
% if ~exist('nSubsample','var')
%     nSubsample = [];
% end

%% inputs 
p = meg_params(sprintf('TANoise_Analysis'));

tIdx = 1:p.trialTime; % epoch time, 1 indexed
t = p.tstart:p.tstop; % epoch time, relative to precue

targetNames = p.eventNames(2:3); 
fields = fieldnames(D4); % {'cueT1','cueT2','cueN'}; % fieldnames(D.subject); 

% peak timing 
switch expt 
    case 'TANoise'
        windowPostTarget = 80; % 70 min window post target to start looking for evoked max, 0 if right after target onset 
        windowSize = 170;  % 190 TANoise(1)  window length to look for target evoked max 
    case 'TA2'
        windowPostTarget = 50; % 100 
        windowSize = 250; % 200 
end
A.windowPostTarget = windowPostTarget;
A.windowSize = windowSize; 
[sessionNames,subjectNames] = meg_sessions(expt); 

%% setup 
plotFigs = 1; 
exptDir = meg_pathToTAMEG(expt, user);
figDir = sprintf('%s/Group/figures/%s_max/nCh_%d_%d',exptDir,slice,nChOI(1),nChOI(end)); 
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% subject scatterplot colors 
colors = distinguishable_colors(numel(subjectNames));

%% session structure: D5.cue(time x session) 
D5 = []; vals = []; 
D5.nChannels = nChOI; 
for iF = 1:numel(fields)
    for iSession = 1:numel(sessionNames)
        channels = groupC(iSession).channelsRanked(nChOI);
        D5.session.channels(iSession,:) = channels; 
        promDir = groupC(iSession).channelDirection(channels);
        D5.session.chPromDir(iSession,:) = promDir; 
        
        %  vals = D4.(fields{iF})(:,nChOI,iSession); % subsampling edit 
        vals = D4.(fields{iF})(:,channels,iSession); 
        
        % vals = vals.*promDir'; % flip data direction by direction of peak prominence
        vals = nanmean(vals,2); % average across channels
        
        D5.session.(fields{iF})(:,:,iSession) = vals;
    end
    D5.session.(fields{iF}) = squeeze(D5.session.(fields{iF})); 
end

%% subject structure: D5.cue(time x subject) 
D5.subject = []; vals = []; 
nSubject = 1; 
for iSession = 1:2:numel(sessionNames)
    for iF = 1:numel(fields)
         vals = nanmean(D5.session.(fields{iF})(:,iSession:iSession+1),2); 
         D5.subject.(fields{iF})(:,nSubject) = vals; 
    end
    nSubject = nSubject + 1; 
end

%% group structure: D5.cue(time x 1)
D5.group = []; vals = [];
for iF = 1:numel(fields)
    D5.group.(fields{iF}) = nanmean(D5.subject.(fields{iF}),2); 
end

%% exchange variable names 
D = D5; 

%% setup A structure 
%  stores max amp and time
%  by A.targetName.cueName: trials x ch x subjects x max
A = []; M = []; I = []; 
for iT = 1:numel(targetNames)
    for iF = 1:numel(fields) 
        A.(targetNames{iT}).(fields{iF}) = [];
    end
end

%% find max and max idx per target per condition per subject 
for iF = 1:numel(fields) % cue condition
    for iT = 1:numel(targetNames) % target
        for iSubject = 1:numel(subjectNames) % subject
            target = targetNames(iT);
            idxTargetName = find(strcmp(p.eventNames,target));
            idxT = find(t==p.eventTimes(idxTargetName));
            window = idxT+windowPostTarget:idxT+windowSize;
            
            % by subject
            [M(iSubject),I(iSubject)] = max(D.subject.(fields{iF})(window,iSubject));
            A.(targetNames{iT}).(fields{iF}).max = M;
            A.(targetNames{iT}).(fields{iF}).maxIdx = t(window(I));
            
        end
        % group avg
        [groupM,groupI] = max(D.group.(fields{iF})(window));
        A.(targetNames{iT}).(fields{iF}).groupMax = groupM;
        A.(targetNames{iT}).(fields{iF}).groupMaxIdx = t(window(groupI));
    end
end

targetNames = {'T1','T2'}; 
for iF = 1:numel(fields) % cue condition
    for iT = 1:numel(targetNames) % target
            target = targetNames(iT);
            idxTargetName = find(strcmp(p.eventNames,target));
            idxT = find(t==p.eventTimes(idxTargetName));
            window = idxT+windowPostTarget:idxT+windowSize;
        for iSession = 1:numel(sessionNames) % subject
            % by session
            [M(iSession),I(iSession)] = max(D.session.(fields{iF})(window,iSession));
            A.(targetNames{iT}).(fields{iF}).maxSession = M;
            A.(targetNames{iT}).(fields{iF}).maxIdxSession = t(window(I));
        end
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

%% make figures!
if plotFigs
    %% fig 1: scatterplot max w v and h error bars
    % find max and max idx per subject, then average
    figure
    set(gcf,'Position',[100 100 1000 400])
    
    for iT = 1:numel(targetNames)
        target = targetNames(iT);
        idxTargetName = find(strcmp(p.eventNames,target));
        idxT = find(t==p.eventTimes(idxTargetName));
        window = idxT:idxT+windowPostTarget+windowSize;
        
        subplot (1,numel(targetNames),iT)
        meg_figureStyle
        hold on
        for iF = 1:numel(fields)
            % vertical error bar
            errV = errorbar(A.(targetNames{iT}).(fields{iF}).maxIdxMean,... % X: max idx
                A.(targetNames{iT}).(fields{iF}).maxMean,... % Y: max
                A.(targetNames{iT}).(fields{iF}).maxStd/numel(subjectNames),... % vertical error bar: max ste
                'Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',0);
            % horizontal error bar
            errV = errorbar(A.(targetNames{iT}).(fields{iF}).maxIdxMean,... % X: max idx
                A.(targetNames{iT}).(fields{iF}).maxMean,... % Y: max
                A.(targetNames{iT}).(fields{iF}).maxIdxStd/numel(subjectNames),... % vertical error bar: max ste
                'horizontal','Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',0);
            
            p1 = plot(A.(targetNames{iT}).(fields{iF}).maxIdxMean,...
                A.(targetNames{iT}).(fields{iF}).maxMean,...
                'o','MarkerSize',5);
            p1.Color = p.cueColors(iF,:);
            p1.MarkerFaceColor = p.cueColors(iF,:);
            p.cueColors(iF,:);
        end
        % ylim([0e-13,2.6e-13])
        xlim([t(window(1))-10,t(window(end))])
        xlabel('Time (ms)')
        ylabel('Max amp (T)')
        vline(p.eventTimes(iT+1),'k:','target')
        title(sprintf('%s',targetNames{iT}))
    end
    sgtitle(sprintf('max by cue condition (channels = %d:%d, selection = %s)',nChOI(1),nChOI(end),groupC(1).selectionMethod))
    
    %% fig 2: scatterplot max w v error bars
    % group max 
    figure
    set(gcf,'Position',[100 100 1000 400])
        
    for iT = 1:numel(targetNames)
        target = targetNames(iT);
        idxTargetName = find(strcmp(p.eventNames,target));
        idxT = find(t==p.eventTimes(idxTargetName));
        window = idxT:idxT+windowPostTarget+windowSize;
        
        subplot (1,numel(targetNames),iT)
        meg_figureStyle
        hold on
        for iF = 1:numel(fields)
            % vertical error bar
            errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
                A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
                A.(targetNames{iT}).(fields{iF}).maxStd/numel(subjectNames),... % vertical error bar: max ste
                'Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',0);
            % horizontal error bar
%             errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
%                 A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
%                 A.(targetNames{iT}).(fields{iF}).maxIdxStd/numel(subjectNames),... % vertical error bar: max ste
%                 'horizontal','Color',p.cueColors(iF,:),'LineWidth',2,'CapSize',0);
            
            p1 = plot(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,...
                A.(targetNames{iT}).(fields{iF}).groupMax,...
                'o','MarkerSize',5);
            p1.Color = p.cueColors(iF,:);
            p1.MarkerFaceColor = p.cueColors(iF,:);
            p.cueColors(iF,:);
        end
        % ylim([0e-13,2.6e-13]) % ylim([0.5e-13,2.6e-13])
        xlim([t(window(1))-10,t(window(end))])
        xlabel('time (ms)')
        ylabel('amp (T)')
        vline(p.eventTimes(iT+1),'k:','target')
        title(sprintf('%s',targetNames{iT}))
    end
    sgtitle(sprintf('max (first average group TS) by cue condition (channels = %d:%d, selection = %s)',nChOI(1),nChOI(end),groupC(1).selectionMethod))
    
    %% fig 3.1: group level plot ERF w max scatter
    % max averaged across subject ERFs 
    p.colorTarget = [0.5 0.5 0.5; 0 0 0];
    
    figure
    meg_figureStyle
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
            window = idxT+windowPostTarget:idxT+windowSize;
            % vertical error bar
            errV = errorbar(A.(targetNames{iT}).(fields{iF}).maxIdxMean,... % X: max idx
                A.(targetNames{iT}).(fields{iF}).maxMean,... % Y: max
                A.(targetNames{iT}).(fields{iF}).maxStd/numel(subjectNames),... % vertical error bar: max ste
                'Color',p.cueColors(iF,:),'LineWidth',1,'CapSize',0);
            % horizontal error bar
            errV = errorbar(A.(targetNames{iT}).(fields{iF}).maxIdxMean,... % X: max idx
                A.(targetNames{iT}).(fields{iF}).maxMean,... % Y: max
                A.(targetNames{iT}).(fields{iF}).maxIdxStd/numel(subjectNames),... % vertical error bar: max ste
                'horizontal','Color',p.cueColors(iF,:),'LineWidth',1,'CapSize',0);
            
            % plot scatter point of max avg by max idx
            p1 = plot(A.(targetNames{iT}).(fields{iF}).maxIdxMean,...
                A.(targetNames{iT}).(fields{iF}).maxMean,...
                'o','MarkerSize',5);
            p1.Color = p.cueColors(iF,:);
            p1.MarkerFaceColor = p.cueColors(iF,:);
            
            % plot time window to find max
            plot(t(window),zeros(size(window)),'LineWidth',4,'Color',p.colorTarget(iT,:))
        end
    end
    yline(0,'k--')
    xlim([t(1),t(end)])
    xlabel('time (ms)')
    ylabel('amp (T)')
    xline(p.eventTimes,'k--',p.eventNames)
    ax1 = gca; 
    sgtitle(sprintf('group max by target by cue condition (nChannels = %d:%d)',nChOI(1),nChOI(end)))
    
    %% fig 3.2: group level plot ERF 
    p.colorTarget = [0.5 0.5 0.5; 0 0 0];
    
    figure
    meg_figureStyle
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
            window = idxT+windowPostTarget:idxT+windowSize;
            % vertical error bar
%             errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
%                 A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
%                 A.(targetNames{iT}).(fields{iF}).maxStd/numel(subjectNames),... % vertical error bar: max ste
%                 'Color',p.cueColors(iF,:),'LineWidth',1,'CapSize',0);
%             % horizontal error bar
%             errV = errorbar(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,... % X: max idx
%                 A.(targetNames{iT}).(fields{iF}).groupMax,... % Y: max
%                 A.(targetNames{iT}).(fields{iF}).maxIdxStd/numel(subjectNames),... % vertical error bar: max ste
%                 'horizontal','Color',p.cueColors(iF,:),'LineWidth',1,'CapSize',0);
            
            % plot scatter point of max avg by max idx
%             p1 = plot(A.(targetNames{iT}).(fields{iF}).groupMaxIdx,...
%                 A.(targetNames{iT}).(fields{iF}).groupMax,...
%                 'o','MarkerSize',5);
%             p1.Color = p.cueColors(iF,:);
%             p1.MarkerFaceColor = p.cueColors(iF,:);
            
            % plot time window to find max
            plot(t(window),zeros(size(window)),'LineWidth',4,'Color',p.colorTarget(iT,:))
        end
    end
    yline(0,'k--')
    xlim([t(1),t(end)])
    ylim(ax1.YLim); 
    xlabel('time (ms)')
    ylabel('amp (T)')
    xline(p.eventTimes,'k--',p.eventNames)
    title(sprintf('group max by target by cue condition (nChannels = %d:%d)',nChOI(1),nChOI(end)))
    
    %% fig 4: subject level plot ERF w max scatter
    p.colorTarget = [0.5 0.5 0.5; 0 0 0];
    
    figure
    meg_figureStyle
    set(gcf,'Position',[100 100 1400 800])
    for iSubject = 1:numel(subjectNames)
        subplot (5,2,iSubject)
        hold on
        for iF = numel(fields):-1:1
            % plot ERF
            p1 = plot(t,D.subject.(fields{iF})(:,iSubject),'LineWidth',1);
            p1.Color = p.cueColors(iF,:);
            p1.Color(4) = p.colorAlpha;
            for iT = 1:numel(targetNames)
                target = targetNames(iT);
                idxTargetName = find(strcmp(p.eventNames,target));
                idxT = find(t==p.eventTimes(idxTargetName));
                window = idxT+windowPostTarget:idxT+windowSize;
                
                % plot scatter point of max by max idx
                p1 = plot(A.(targetNames{iT}).(fields{iF}).maxIdx(iSubject),...
                    A.(targetNames{iT}).(fields{iF}).max(iSubject),...
                    'o','MarkerSize',5);
                p1.Color = p.cueColors(iF,:);
                p1.MarkerFaceColor = p.cueColors(iF,:);
                
                % plot time window to find max
                plot(t(window),zeros(size(window)),'LineWidth',4,'Color',p.colorTarget(iT,:))
            end
        end
        yline(0,'k--')
        xlim([t(1),t(end)])
        xlabel('time (ms)')
        ylabel('amp (T)')
        xline(p.eventTimes,'k--',p.eventNames)
        title(sprintf('%s',subjectNames{iSubject}))
    end
    % rd_supertitle2(sprintf('(nChannels = %d)',nChOI))
    
    %% fig 5: session level plot ERF w max scatter
    
    p.colorTarget = [0.5 0.5 0.5; 0 0 0];
    
    figure
    set(gcf,'Position',[100 100 600 1200])
    for iSession = 1:numel(sessionNames)
        subplot (numel(subjectNames),2,iSession)
        hold on
        for iF = numel(fields):-1:1
            % plot ERF
            p1 = plot(t,D.session.(fields{iF})(:,iSession),'LineWidth',1);
            p1.Color = p.cueColors(iF,:);
            p1.Color(4) = p.colorAlpha;
            for iT = 1:numel(targetNames)
                target = targetNames(iT);
                idxTargetName = find(strcmp(p.eventNames,target));
                idxT = find(t==p.eventTimes(idxTargetName));
                window = idxT+windowPostTarget:idxT+windowSize;
                %
                %             % plot scatter point of max by max idx
                %             p1 = plot(A.(targetNames{iT}).(fields{iF}).maxIdx(iSubject),...
                %                 A.(targetNames{iT}).(fields{iF}).max(iSubject),...
                %                 'o','MarkerSize',5);
                %             p1.Color = p.cueColors(iF,:);
                %             p1.MarkerFaceColor = p.cueColors(iF,:);
                
                % plot time window to find max
                plot(t(window),zeros(size(window)),'LineWidth',4,'Color',p.colorTarget(iT,:))
            end
        end
        yline(0,'k--')
        xlim([t(1),t(end)])
        xline(p.eventTimes,'k--')
        title(sprintf('%s',und2space(sessionNames{iSession})))
    end
    xlabel('time (ms)')
    ylabel('amp (T)')
    % rd_supertitle2(sprintf('(nChannels = %d)',nChOI))
    
    %% fig 6: subject scatter compare max by cond by target
    
    % fig 6.1: cue T1 v cue T2 
    
    figure
    set(gcf,'Position',[0 0 1000 400])
    for iT = 1:numel(targetNames)
        subplot (1,numel(targetNames),iT)
        hold on
        for iSubject = 1:numel(subjectNames)
            conds = {'cueT1','cueT2'}; % 'cueN'
            p1 = plot(A.(targetNames{iT}).(conds{1}).max(iSubject),... % cue T1
                A.(targetNames{iT}).(conds{2}).max(iSubject),... % cue T2
                'o','MarkerSize',15);
            p1.Color = 'k';
            p1.MarkerFaceColor = colors(iSubject,:);
        end
        axis equal
        l = refline(1,0);
        l.Color = 'k';
        l.LineStyle = '--';
        xlabel(sprintf('%s',conds{1}))
        ylabel(sprintf('%s',conds{2}))
        title(sprintf('%s max (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
    end
    
%     % fig 6.2: cue T1 v cue N 
%     figure
%     set(gcf,'Position',[0 0 1000 400])
%     for iT = 1:numel(targetNames)
%         subplot (1,numel(targetNames),iT)
%         hold on
%         for iSubject = 1:numel(subjectNames)
%             conds = {'cueT1','cueN'}; % 'cueN'
%             p1 = plot(A.(targetNames{iT}).(conds{1}).max(iSubject),... % cue T1
%                 A.(targetNames{iT}).(conds{2}).max(iSubject),... % cue T2
%                 'o','MarkerSize',15);
%             p1.Color = 'k';
%             p1.MarkerFaceColor = colors(iSubject,:);
%         end
%         axis equal
%         l = refline(1,0);
%         l.Color = 'k';
%         l.LineStyle = '--';
%         xlabel(sprintf('%s',conds{1}))
%         ylabel(sprintf('%s',conds{2}))
%         title(sprintf('%s max (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
%     end
%     
%     % fig 6.3: cue T2 v cue N 
%     figure
%     set(gcf,'Position',[0 0 1000 400])
%     for iT = 1:numel(targetNames)
%         subplot (1,numel(targetNames),iT)
%         hold on
%         for iSubject = 1:numel(subjectNames)
%             conds = {'cueT2','cueN'}; % 'cueN'
%             p1 = plot(A.(targetNames{iT}).(conds{1}).max(iSubject),... % cue T1
%                 A.(targetNames{iT}).(conds{2}).max(iSubject),... % cue T2
%                 'o','MarkerSize',15);
%             p1.Color = 'k';
%             p1.MarkerFaceColor = colors(iSubject,:);
%         end
%         axis equal
%         l = refline(1,0);
%         l.Color = 'k';
%         l.LineStyle = '--';
%         xlabel(sprintf('%s',conds{1}))
%         ylabel(sprintf('%s',conds{2}))
%         title(sprintf('%s max (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
%     end
    
    %% fig 7: subject scatter compare max idx by cond by target
    
    % fig 7.1: cue T1 v cue T2 
    figure
    set(gcf,'Position',[0 0 1000 400])
    for iT = 1:numel(targetNames)
        subplot (1,numel(targetNames),iT)
        hold on
        target = targetNames(iT);
        idxTargetName = find(strcmp(p.eventNames,target));
        idxT = find(t==p.eventTimes(idxTargetName));
        window = idxT:idxT+windowPostTarget+windowSize;
        for iSubject = 1:numel(subjectNames)
            conds = {'cueT1','cueT2'};
            p1 = plot(A.(targetNames{iT}).(conds{1}).maxIdx(iSubject),... % cue T1
                A.(targetNames{iT}).(conds{2}).maxIdx(iSubject),... % cue T2
                'o','MarkerSize',15);
            p1.Color = 'k';
            p1.MarkerFaceColor = colors(iSubject,:);
        end
        xlim([t(window(1)),t(window(end))])
        ylim([t(window(1)),t(window(end))])
        l = refline(1,0);
        l.Color = 'k';
        l.LineStyle = '--';
        xlabel(sprintf('%s',conds{1}))
        ylabel(sprintf('%s',conds{2}))
        title(sprintf('%s max idx (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
    end
    
%     % fig 7.2: cue T1 v cue N 
%     figure
%     set(gcf,'Position',[0 0 1000 400])
%     for iT = 1:numel(targetNames)
%         subplot (1,numel(targetNames),iT)
%         hold on
%         target = targetNames(iT);
%         idxTargetName = find(strcmp(p.eventNames,target));
%         idxT = find(t==p.eventTimes(idxTargetName));
%         window = idxT:idxT+windowPostTarget+windowSize;
%         for iSubject = 1:numel(subjectNames)
%             conds = {'cueT1','cueN'};
%             p1 = plot(A.(targetNames{iT}).(conds{1}).maxIdx(iSubject),... % cue T1
%                 A.(targetNames{iT}).(conds{2}).maxIdx(iSubject),... % cue T2
%                 'o','MarkerSize',15);
%             p1.Color = 'k';
%             p1.MarkerFaceColor = colors(iSubject,:);
%         end
%         xlim([t(window(1)),t(window(end))])
%         ylim([t(window(1)),t(window(end))])
%         l = refline(1,0);
%         l.Color = 'k';
%         l.LineStyle = '--';
%         xlabel(sprintf('%s',conds{1}))
%         ylabel(sprintf('%s',conds{2}))
%         title(sprintf('%s max idx (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
%     end
%     
%     % fig 7.3: cue T2 v cue N 
%     figure
%     set(gcf,'Position',[0 0 1000 400])
%     for iT = 1:numel(targetNames)
%         subplot (1,numel(targetNames),iT)
%         hold on
%         target = targetNames(iT);
%         idxTargetName = find(strcmp(p.eventNames,target));
%         idxT = find(t==p.eventTimes(idxTargetName));
%         window = idxT:idxT+windowPostTarget+windowSize;
%         for iSubject = 1:numel(subjectNames)
%             conds = {'cueT2','cueN'};
%             p1 = plot(A.(targetNames{iT}).(conds{1}).maxIdx(iSubject),... % cue T1
%                 A.(targetNames{iT}).(conds{2}).maxIdx(iSubject),... % cue T2
%                 'o','MarkerSize',15);
%             p1.Color = 'k';
%             p1.MarkerFaceColor = colors(iSubject,:);
%         end
%         xlim([t(window(1)),t(window(end))])
%         ylim([t(window(1)),t(window(end))])
%         l = refline(1,0);
%         l.Color = 'k';
%         l.LineStyle = '--';
%         xlabel(sprintf('%s',conds{1}))
%         ylabel(sprintf('%s',conds{2}))
%         title(sprintf('%s max idx (nChannels = %d:%d)',targetNames{iT},nChOI(1),nChOI(end)))
%    end
    
    %% topo channel(s) selected and prom (fix not working w new fieldtrip) 
    %% session 1
%     figure
%     set(gcf,'Position',[0 0 1600 500])
%     nFig = 1; 
%     for iS = 1:2:numel(sessionNames)
%         subplot(1,numel(subjectNames),nFig)
%         hold on 
%         chs = []; vals = []; 
%         
%         chs = groupC(iS).channelsRanked(nChOI); 
%         promDir = groupC(iS).channelDirection(chs); 
%         
%         vals = zeros(157,1); 
%         vals(chs) = promDir; 
%         meg_topoplot(vals)
%         caxis([-1,1])
%         % title(sprintf('%s',und2space(sessionNames{iS})))
%         nFig = nFig + 1; 
%     end
% 
%     %% session 2 
%     figure
%     set(gcf,'Position',[0 0 1600 500])
%     nFig = 1; 
%     for iS = 2:2:numel(sessionNames)
%         subplot(1,numel(subjectNames),nFig)
%         chs = []; vals = []; 
%         
%         chs = groupC(iS).channelsRanked(nChOI); 
%         promDir = groupC(iS).channelDirection(chs); 
%         
%         vals = zeros(157,1); 
%         vals(chs) = promDir; 
%         meg_topoplot(vals)
%         caxis([-1,1])
%         % title(sprintf('%s',und2space(sessionNames{iS})))
%         nFig = nFig + 1; 
%     end
    
    %% save figs
    fH = sort(double(findobj(0,'Type','figure')));
    rd_saveAllFigs(fH,...
        {'targetMax_subject','targetMax_group','groupERF_subject','groupERF_group','subjectERF','sessionERF',...
        'scatter_subject_cT1cT2_max',... % 'scatter_subject_cT1cN_max','scatter_subject_cT2cN_max',...
        'scatter_subject_cT1cT2_maxIdx'},... %'scatter_subject_cT1cN_maxIdx','scatter_subject_cT2cN_maxIdx'...        % 'topo_session1','topo_session2'},...
        sprintf('ch%d_%d',nChOI(1),nChOI(end)),figDir);
    
    %% save analysis mat 
    D = D5; 
    save(sprintf('%s/ERFbyCue.mat',figDir),'D')
    save(sprintf('%s/max.mat',figDir),'A')
  
end


