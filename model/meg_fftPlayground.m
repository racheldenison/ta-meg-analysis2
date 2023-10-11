function meg_fftPlayground

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian (personal), karen (lab)

addpath(genpath(pwd))

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPCfit_separate1/%s',user,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Prepare data (separate precue T1 and T2) 
% --- Load ITPC data ---
load(sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/unused/groupA_ITPCspectrogram_byAtt.mat',user)) % A variable 
% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

% --- Data settings ---
% sampling = 1:1:7001; % 10 ms 
foi = 20; % frequency of interest, Hz
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 

Fs = 1000; % sampling frequency 

fitLevel = 'session'; 
baselineToi = abs(p.tstart)+p.eventTimes(1)-100:abs(p.tstart)+p.eventTimes(1); 

% -- MEG settings --- 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

% --- Extract data into precue conds at desired foi --- 
cueNames = {'all','cueT1','cueT2'};
fitLevels = {'session','subject','group'}; 
clear data 
for i = 1:numel(cueNames)
    data.(cueNames{i}).session = squeeze(A.(cueNames{i}).session(foi,:,:)); 
    data.(cueNames{i}).subject = squeeze(A.(cueNames{i}).subject(foi,:,:)); % frequency x time x subjects
    data.(cueNames{i}).group   = squeeze(mean(data.(cueNames{i}).subject,2)); % average across subjects
end

%% Concatenate all data into t x 40 (sessions and precue) 
dataAll = cat(2,data.(cueNames{1}).session,data.(cueNames{2}).session); 
% extract frequency of interest 
dataAll = dataAll(tIdx,:);
nSessions = size(dataAll,2);

%% FFT on data 
user = 'kantian';
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPCfit_separate1/FFT/%s_unpadded',user,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

clear f amp
for i = 1:nSessions
    [f(i,:),amp(i,:)] = meg_fft(dataAll(:,i));
    sgtitle(sprintf('FFT %d',i))
    figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_%d',i);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Plot group averaged fft
meanF = mean(f,1,'omitnan'); 
meanAmp = mean(amp,1,'omitnan'); 

figure
subplot 121
hold on
meg_figureStyle
plot(meanF,meanAmp,'LineWidth',2,'Color',p.cueColors(1,:))
ylabel('Amplitude')
xlabel('Frequency (Hz)')
xlim([0 5])

subplot 122
hold on
meg_figureStyle
plot(log(meanF),log(meanAmp),'LineWidth',2,'Color',p.cueColors(1,:))
ylabel('Log(Amplitude)')
xlabel('Log(Frequency) (Hz)')
xlim([0 5])

sgtitle(sprintf('FFT Group n = %d',nSessions))
figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_group_n%d',nSessions);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% FFT on 2 Hz sine (simulated) 
freq = 2; 
t = 1:2048;
test = sin( (freq*pi/(Fs/2)) * t); 
[fTest,ampTest] = meg_fft(test');
sgtitle(sprintf('FFT simulation'))
figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_2048');
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% Bootstrap ITPC FFTs 
nBoot = 1000; 
for iB = 1:nBoot
    idx = randi(40,[1 nSessions]);
    for iS = 1:nSessions
        dataBoot(iS,:,iB) = amp(idx(iS),:); 
    end
end

%% Plot bootstrapped group FFT 
figure
set(gcf,'Position',[100 100 800 400])
for i = 1:nBoot
    val = mean(dataBoot(:,:,i),1,'omitnan'); 
    subplot 121
    hold on
    meg_figureStyle
    plot(f,val,'LineWidth',0.5,'Color',p.cueColors(1,:)) % bootstrapped data 
    plot(meanF,meanAmp,'LineWidth',2,'Color','k') % real data 
    ylabel('Amplitude')
    xlabel('Frequency (Hz)')
    xlim([0 5])
    xline(2)

    subplot 122
    hold on
    meg_figureStyle
    plot(log(f),log(val),'LineWidth',0.5,'Color',p.cueColors(1,:)) % bootstrapped data 
    plot(log(meanF),log(meanAmp),'LineWidth',2,'Color','k') % real data 
    ylabel('Log(Amplitude)')
    xlabel('Log(Frequency) (Hz)')
    xlim([0 5])
    xline(log(2))
end

% fit line to start and end? 
x = [log(meanF(2)) log(meanF(end))];
y = [log(meanAmp(2)) log(meanAmp(end))];
plot(x,y,'r')

sgtitle(sprintf('FFT Bootstrapped (n = %d, nBoots = %d)',nSessions,nBoot))
figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_group_n%d_bootstrapped',nSessions);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))




