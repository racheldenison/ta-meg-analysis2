% meg_plotSubjectBehav

% MEG_PLOTSUBJECTBEHAV
%
% Karen Tian
% August 2020

%% input 
expt = 'TANoise';
saveFigs = 1; 
figDir = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/ITPCUpDown'; 

%% setup
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); % session & subject info 
p = meg_params(sprintf('%s_Analysis',expt)); % expt parameters

targets = {'T1','T2'};
cue = {'cueT1','cueT2'}; 
ITPCDirs = [1,-1]; % up, down
colors = [187/255, 213/255, 151/255;... 
    190/255, 151/255, 198/255]; % ITPC up, down 
dotSize= 180; 

%% load behavior 
load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/mat_byCue/groupB.mat') % groupB

%% plot dprime by uppers downers
figure
hold on 
for iT = 1:numel(targets)
    for iC = 1:numel(cue)
        if iT == iC
            spacer = 1; % valid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectDprime;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        else
            spacer = 2; % invalid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectDprime;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        end
        allyVals(iT,iC,:) = yVals;
        allxVals(iT,iC,:) = xVals;
    end
end

% visualize cueing effect, connect subjects 
for iT = 1:numel(targets)
    for iS = 1:numel(yVals)
        if ITPCsubject(iS) ~= 0
            plot([allxVals(iT,1,iS),allxVals(iT,2,iS)],...
                [allyVals(iT,1,iS),allyVals(iT,2,iS)],...
                'Color',[0.5 0.5 0.5 0.5],'LineWidth',1)
        end
    end
end

for iT = 1:numel(targets)
    for iC = 1:numel(cue)
        if iT == iC
            spacer = 1; % valid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectDprime;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        else
            spacer = 2; % invalid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectDprime;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        end
        for iS = 1:numel(yVals)
            if ITPCsubject(iS) == 1 % upper
                scatter(xVals(iS),yVals(iS),dotSize,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1)
            elseif ITPCsubject(iS) == -1 % downer
                scatter(xVals(iS),yVals(iS),dotSize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1],'LineWidth',1)
            end
        end
    end
end

xlim([0.5,3.5])
ylabel('Sensitivity (d'')')
ax = gca; 
ax.FontSize = 12; 
ax.LineWidth = 1.5; 
ax.TickDir = 'out'; 

if saveFigs 
    titleText = 'dprime'; 
    print(sprintf('%s/behavior/%s.eps',figDir,titleText),'-depsc')
end

%% plot rt by uppers downers
figure
hold on 

for iT = 1:numel(targets)
    for iC = 1:numel(cue)
        if iT == iC
            spacer = 1; % valid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectRT;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        else
            spacer = 2; % invalid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectRT;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        end
        allyVals(iT,iC,:) = yVals;
        allxVals(iT,iC,:) = xVals;
    end
end

% visualize cueing effect, connect subjects 
for iT = 1:numel(targets)
    for iS = 1:numel(yVals)
        if ITPCsubject(iS) ~= 0
            plot([allxVals(iT,1,iS),allxVals(iT,2,iS)],...
                [allyVals(iT,1,iS),allyVals(iT,2,iS)],...
                'Color',[0.5 0.5 0.5 0.5],'LineWidth',1)
        end
    end
end

for iT = 1:numel(targets)
    for iC = 1:numel(cue)
        if iT == iC
            spacer = 1; % valid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectRT;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        else
            spacer = 2; % invalid
            yVals = groupB.(targets{iT}).(cue{iC}).subjectRT;
            xVals = zeros(1,size(yVals,2))+iT+spacer/5;
        end
        for iS = 1:numel(yVals)
            if ITPCsubject(iS) == 1 % upper
                scatter(xVals(iS),yVals(iS),dotSize,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',1)
            elseif ITPCsubject(iS) == -1 % downer
                scatter(xVals(iS),yVals(iS),dotSize,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1],'LineWidth',1)
            end
        end
    end
end

xlim([0.5,3.3])
ylabel('Reaction time (ms)')
ax = gca; 
ax.FontSize = 12; 
ax.LineWidth = 1.5; 
ax.TickDir = 'out'; 

if saveFigs 
    titleText = 'rt'; 
    print(sprintf('%s/behavior/%s.eps',figDir,titleText),'-depsc')
end

%% cueT1 cueT2 scatter
for iT = 1:numel(targets)
    figure
    hold on
    xVals = groupB.(targets{iT}).cueT1.subjectDprime; % cueT1
    yVals = groupB.(targets{iT}).cueT2.subjectDprime; % cueT2
    for iS = 1:numel(subjectNames)
        if ITPCsubject(iS) == 1 % upper
            scatter(xVals(iS),yVals(iS),dotSize,'filled','MarkerFaceColor',colors(1,:))
        elseif ITPCsubject(iS) == -1 % downer
            scatter(xVals(iS),yVals(iS),dotSize,'filled','MarkerFaceColor',colors(2,:))
        end
    end
    xlim([0 3.5])
    ylim([0 3.5])
    ref = refline(1,0);
    ref.Color = [0.7 0.7 0.7]; 
    ref.LineStyle = '--'; 
    title(sprintf('%s Sensitivity',targets{iT}))
    xlabel('Cue T1')
    ylabel('Cue T2')
    axis square
    ax = gca; 
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    if saveFigs
        titleText = sprintf('%s',targets{iT});
        print(sprintf('%s/behavior/%s.eps',figDir,titleText),'-depsc')
    end
end




