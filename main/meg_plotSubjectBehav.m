% meg_plotSubjectBehav

% MEG_PLOTSUBJECTBEHAV
%
% Karen Tian
% August 2020

%% input 
expt = 'TANoise';
saveFigs = 1; 

%% setup
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); % session & subject info 
p = meg_params(sprintf('%s_Analysis',expt)); % expt parameters

targets = {'T1','T2'};
cue = {'cueT1','cueT2'}; 
ITPCDirs = [1,-1]; % up, down
colors = [97/255, 159/255, 145/255;... 
    179/255, 117/255, 195/255]; % ITPC up, down 

%% load behavior 
load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/mat_byCue/groupB.mat') % groupB

%% plot dprime by uppers downers
figure
hold on 
for iT = 1:numel(targets)
    for iC = 1:numel(cue)
        yVals = groupB.(targets{iT}).(cue{iC}).subjectDprime; 
        xVals = zeros(1,size(yVals,2))+iT+iC/5; 
        for iS = 1:numel(yVals)
            if ITPCsubject(iS) == 1 % upper
                scatter(xVals(iS),yVals(iS),50,'filled','MarkerFaceColor',colors(1,:))
            elseif ITPCsubject(iS) == -1 % downer
                scatter(xVals(iS),yVals(iS),50,'filled','MarkerFaceColor',colors(2,:))
            end
        end
    end
end
xlim([0,3.5])

if saveFigs 
    titleText = 'phaseAngle_T1'; 
    print(sprintf('%s/%s.eps',figDir,titleText),'-depsc')
end
