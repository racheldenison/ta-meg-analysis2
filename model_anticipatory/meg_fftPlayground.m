function meg_fftPlayground

%% Settings
saveFigs = 0; 
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
    if saveFigs 
        saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    end
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

if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

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
        dataBoot(iS,:,iB) = amp(idx(iS),:); % sessions/precue  x freqs x boots 
    end
end

% Average bootstrapped data (along sessions/precue) to group
dataBootGroup = squeeze(mean(dataBoot,1,'omitnan')); % freqs x boots 

% Calculate 5th and 95th percentile 
for iAmp = 1:size(dataBootGroup,1) % loop amps 
    low = 5/2; % 2.5 percentile 
    high = 100-5/2; % 97.5 percentile 
    pLowBoot(iAmp) = prctile(dataBootGroup(iAmp,:),low);
    pHighBoot(iAmp) = prctile(dataBootGroup(iAmp,:),high);
    pLowBootLog(iAmp) = prctile(log(dataBootGroup(iAmp,:)),low);
    pHighBootLog(iAmp) = prctile(log(dataBootGroup(iAmp,:)),high);
end

f = f(1,:); % extract just 1st row of frequencies from fft (should be same for all 40 session/precue) 

%% Fit 1/f to log-log freq amp 
% Linear fit from 35 to 200 Hz
% See Hermes, Miller, Wandell & Winawer 2014 
fStart = 35; 
fEnd = 200; 
fSkip = 60:60:fEnd; 
[~,fStart_idx] = min( abs( fStart - meanF ) );
[~,fEnd_idx] = min( abs( fEnd - meanF ) );
for i = 1:numel(fSkip)
    [~,fSkip_idx(i)] = min( abs( fSkip(i) - meanF ) );
end
fSkip_idx = fSkip_idx - fStart;

freqToFit = log( meanF(fStart_idx:fEnd_idx   ));
ampToFit  = log( meanAmp(fStart_idx:fEnd_idx ));

% Remove line noise and harmonics
freqToFit(fSkip_idx) = NaN; 
ampToFit(fSkip_idx) = NaN; 

% Linear fit 
fitP = polyfit( freqToFit(~isnan(freqToFit)), ampToFit(~isnan(ampToFit)) ,1); % slope, intercept 

%% Plot bootstrapped group FFT 
figure
set(gcf,'Position',[100 100 600 300])
bootColor = [0.8 0.8 0.8]; 

% --- Plot f x amp --- 
subplot 121
hold on
meg_figureStyle
plot(f,mean(dataBootGroup,2,'omitnan'),'LineWidth',0.5,'Color',bootColor) % bootstrapped data mean
plot(f,pLowBoot,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 5th percentile 
plot(f,pHighBoot,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 95th percentile 
patch([f fliplr(f)], [pLowBoot fliplr(pHighBoot)],bootColor ,'EdgeColor',bootColor) % shade bootstrapped 5-95th percentile 
plot(meanF,meanAmp,'LineWidth',1,'Color','k') % real data
ylabel('Amplitude')
xlabel('Frequency (Hz)')
xlim([0 5])
% Plot vertical lines at returned integers
for i = 2:11
    xline(f(i))
end
ax = gca;
set(ax ,'Layer', 'Top')

% --- Plot log f x log amp --- 
subplot 122
hold on
meg_figureStyle
logf = log(f); 
plot(logf,log(mean(dataBootGroup,2,'omitnan')),'LineWidth',0.5,'Color',bootColor) % bootstrapped data mean
plot(logf,pLowBootLog,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 5th percentile 
plot(logf,pHighBootLog,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 95th percentile 
patch([logf(2:end) fliplr(logf(2:end))], [pLowBootLog(2:end) fliplr(pHighBootLog(2:end))],bootColor,'EdgeColor',bootColor) % shade bootstrapped 5-95th percentile 
plot(log(meanF),log(meanAmp),'LineWidth',1,'Color','k') % real data
ylabel('Log (amplitude)')
xlabel('Log (frequency)')
xlim([0 log(200)])
% Plot vertical lines at returned integers
for i = 2:11
    xline(logf(i))
end
xline(logf(20)) % check 20 Hz 
xline(logf(60)) % check line noise? 

% Linear fit to log-log 1/f 
x = [log(meanF(fStart_idx)) log(meanF(fEnd_idx))];
y = [log(meanAmp(fStart_idx)) log(meanAmp(fEnd_idx))];
% slope = (y(2)-y(1))/(x(2)-x(1));
% intercept = y(1)-slope*x(1); 
% rLine = refline(slope,intercept); 
x = log(meanF(fStart_idx:fEnd_idx));
y = fitP(1)*x + fitP(2); 
rLine = refline(fitP(1),fitP(2));
rLine.Color = p.cueColors(1,:);
rLine.LineStyle = '-';
plot(x,y,'b','LineWidth',1.5)

% get predicted log amps from log f
fOOF(1) = NaN; 
for i = 2:size(logf,2)
    fOOF(i) = slope*logf(i)+intercept;
end

pOOF(1) = NaN; 
for i = 2:size(f,2)
    pOOF(i) = invprctile( log(dataBootGroup(i,:)) , fOOF(i) );
end

sgtitle(sprintf('FFT Bootstrapped (n = %d, nBoots = %d)',nSessions,nBoot))
figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_group_n%d_bootstrapped',nSessions);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))




