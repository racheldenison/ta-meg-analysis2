% model comparison lrt simulation 

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
amplitude1 = 0.1; 
phase1 = 0; % change here 
% --- Data 2 params --- 
intercept2 = 0.3;
slope2 = 1; % 0.0001; 
amplitude2 = 0.1; 
phase2 = 1.5;
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

%% Plot data 
figure
switch fitLevel
    case 'session'
        set(gcf,'Position',[100 100 500 300])
    case 'subject'
        set(gcf,'Position',[100 100 500 300])
    otherwise
        error('Specify session- or subject-level fit')
end
hold on 

% --- Plot data ---
plot(p.t,dummyData.cueT1,'LineWidth',2)
plot(p.t,dummyData.cueT2,'LineWidth',2)
xlabel('Time (ms)')
ylabel('Simulated ITPC')

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
titleText = sprintf('Separate precue T1 and T2 fits\nintercept1 = %0.2f, slope1 = %0.5f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.5f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
    intercept1,slope1,...
    amplitude1,phase1,...
    intercept2,slope2,...
    amplitude2,phase2,...
    freq);
title(titleText)

ax = gca;
ax.TitleFontSizeMultiplier = 0.7;
ax.TitleFontWeight = 'normal';

% --- Save fig ---
figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_Separate_%s',dateStr);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% Do model fit 









