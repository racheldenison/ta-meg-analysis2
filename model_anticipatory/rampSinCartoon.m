% Phase offsets in two oscillatory signals may manifest as 
% apparent slope differences when assessed by a linear model 
% that does not account for periodicity

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian (persona), karen (lab)

addpath(genpath(pwd))

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPCsimulation_1/%s',user,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

% --- Data settings ---
foi = 20; % frequency of interest, Hz
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 
Fs = 1000; % sampling frequency 
fitLevel = 'session'; 

% -- MEG settings --- 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

% --- Extract data into precue conds at desired foi --- 
cueNames = {'all','cueT1','cueT2'};
fitLevels = {'session','subject','group'}; 

%% --- Generate simulated data ---
Fs = 1000; 
freq = 2; 
% --- Data 1 params --- 
intercept1 = 0.3; 
slope1 = 1; % 0.0001; % 0.0001 more reflects actual data 
amplitude1 = 0.3; 
phase1 = 0; % change here 
% --- Data 2 params --- 
intercept2 = 0.3;
slope2 = 1; % 0.0001; 
amplitude2 = 0.3; 
phase2 = rad2t(pi);
% --- Simulate data (in trial-relative time) --- 
clear dummyData
dummyData.cueT1 = ((slope1/1000 * p.t) + intercept1)  + ( amplitude1 * sin( (freq*pi/(Fs/2)) * (p.t + phase1 * 100 )) ); 
dummyData.cueT2 = ((slope2/1000 * p.t) + intercept2)  + ( amplitude2 * sin( (freq*pi/(Fs/2)) * (p.t + phase2 * 100 )) );

% --- Add noise ---
addNoise = 0; 
if addNoise
    noise1 = (rand(size(dummyData.cueT1))-0.5) * 0.5;
    noise2 = (rand(size(dummyData.cueT2))-0.5) * 0.2;

    dummyData.cueT1 = dummyData.cueT1 + noise1;
    dummyData.cueT2 = dummyData.cueT2 + noise2;
end
%%
figure
set(gcf,'Position',[100 100 500 200])

iF = 1; 
hold on

for iC = 1:numel(cueLevel)
    % --- Data (by precue) ---
    dataFit = dummyData.(cueLevel{iC});

    % --- Plot data ---
    plot(p.t,dataFit,'LineWidth',1,'Color',p.cueColors(iC,:))

    % --- Plot model fit ---
    [minVal,idx] = min(  mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:)  );
    fittedX = mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,:);
    switch fitTypes{iF}
        case 'linear2Hz'
            paramNames = {'intercept','slope','amplitude','phase'};
        case 'linear'
            paramNames = {'intercept','slope'};
    end
    fitColors = {'b','r'};
    plot(t, mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(idx,:),'--','LineWidth',1,'Color',fitColors{iC})
end

% --- Plot event lines ---
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end

% --- Format ---
meg_figureStyle
xlim([-100 2400])
xlabel('Time (ms)')
ylabel('ITPC')

% --- Titles --
switch fitTypes{iF}
    case 'linear2Hz'
        titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(2),...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(3),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(4),...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(2),...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(3),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(4),...
            2);
    case 'linear'
        titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f\nintercept2 = %0.2f, slope2 = %0.3f',...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(2),...
            mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(2));
end
titleText1 = sprintf('Separate precue T1 and T2 fits (%s), fval = %0.2e',...
    mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
    mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval);
titleText_True = sprintf('True: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
    intercept1, slope1, amplitude1, phase1,...
    intercept2, slope2, amplitude2, phase2,...
    freq);
titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted,titleText_True);
title(titleText)

ax = gca;
ax.TitleFontSizeMultiplier = 0.6;
ax.TitleFontWeight = 'normal';

lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast');
lgd.FontSize = 10;

% --- Save fig ---
dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_Separate_%s_%s',fitTypes{iF},dateStr);
saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
