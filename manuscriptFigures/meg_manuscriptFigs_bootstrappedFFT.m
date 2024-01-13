function meg_manuscriptFigs_bootstrappedFFT
% function meg_manuscriptFigs_bootstrappedFFT
% see meg_fftPlayground for full figs

user = 'kantian'; % kantian (personal), karen (lab)

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
annotateStats = 1; 
saveFigs = 0; 
showLegend = 1; 
[style, colors] = meg_manuscriptStyle;
varOI = 'pow'; % amp % plots power or amplitude spectrum 

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('%s/manuscriptFigures/figs',pwd); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Prepare data (separate precue T1 and T2) 
% --- Load ITPC data ---
filename = sprintf('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/groupA_ITPCspectrogram_byAtt.mat',user); 
load(filename) % A variable 
% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

% --- Data settings ---
% sampling = 1:1:7001; % 10 ms 
foi = 20; % frequency of interest, Hz
paddingBefore = 51; % 80 ms before T1 
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

%% Concatenate all data into t x 40 (sessions and precue) or t x 20 (sessions) 
analLevel = 'all'; % 'all' 'precue' 
switch analLevel 
    case 'precue'
        dataAll = cat(2,data.cueT1.session,data.cueT2.session); % oh fuck was this all trials and cue T1
    case 'all'
        dataAll = data.all.session; 
end
% extract frequency of interest
dataAll = dataAll(tIdx,:);
nSessions = size(dataAll,2);

%% FFT on data 
clear f amp pow
for i = 1:nSessions
    [f(i,:),amp(i,:)] = meg_fft(dataAll(:,i));
    % sgtitle(sprintf('FFT %d',i))
    % figTitle = sprintf('TANoise_ITPCFit_Separate_FFT_%d',i);
    % if saveFigs 
    %     saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    % end
end
pow = amp.^2; 

%% Plot group averaged fft
meanF = mean(f,1,'omitnan'); % average sessions 
meanAmp = mean(amp,1,'omitnan'); 
meanPow = mean(pow,1,'omitnan'); 

%% Bootstrap ITPC FFTs 
switch varOI 
    case 'pow'
        var = pow; 
    case 'amp'
        var = amp;
end

nBoot = 1000;
for iB = 1:nBoot
    idx = randi(nSessions,[1 nSessions]);
    for iS = 1:nSessions
        dataBoot(iS,:,iB) = var(idx(iS),:); % sessions/precue  x freqs x boots
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

switch varOI 
    case 'pow'
        varToFit = meanPow; 
    case 'amp'
        varToFit = meanAmp; 
end

freqToFit = log( meanF(fStart_idx:fEnd_idx   ));
varToFit  = log( varToFit(fStart_idx:fEnd_idx ));

% Remove line noise and harmonics
freqToFit(fSkip_idx) = NaN; 
varToFit(fSkip_idx) = NaN; 

% Linear fit 
fitP = polyfit( freqToFit(~isnan(freqToFit)), varToFit(~isnan(varToFit)) ,1); % slope, intercept 

%% Plot bootstrapped group FFT 
figure
set(gcf,'Position',[100 100 312 300])
bootColor = [0.8 0.8 0.8]; 

% --- Plot log f x log amp --- 
fh = subplot (1,1,1); 
hold on
meg_figureStyle
logf = log(f); 
plot(logf,log(mean(dataBootGroup,2,'omitnan')),'LineWidth',0.5,'Color',bootColor) % bootstrapped data mean
plot(logf,pLowBootLog,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 5th percentile 
plot(logf,pHighBootLog,'LineWidth',0.3,'Color',bootColor) % bootstrapped data 95th percentile
pl(1) = patch([logf(2:end) fliplr(logf(2:end))], [pLowBootLog(2:end) fliplr(pHighBootLog(2:end))],bootColor,'EdgeColor',bootColor,'DisplayName','95% CI bootstrapped FFT'); % shade bootstrapped 5-95th percentile
switch varOI 
    case 'pow'
        pl(2) = plot(log(meanF),log(meanPow),'LineWidth',1,'Color','k','DisplayName','Data FFT'); % real data
        ylabel('Log (power)')
        y = [log(meanPow(fStart_idx)) log(meanPow(fEnd_idx))];
    case 'amp'
        pl(2) = plot(log(meanF),log(meanAmp),'LineWidth',1,'Color','k','DisplayName','Data FFT'); % real data
        ylabel('Log (amplitude)')
        y = [log(meanAmp(fStart_idx)) log(meanAmp(fEnd_idx))];
end
xlabel('Frequency (Hz)')
xlim([0 log(200)])
xl = xlim;
yl = ylim;
% Plot vertical lines at returned integers
for i = 2:7
    xh = xline(logf(i),'Color',colors.eventLines); 
    tick(i) = logf(i);
    meg_sendToBack(xh)
    if annotateStats
        txt = meg_annotateStats(tick(i),max(fh.YLim)*1.05,'*'); 
    end
end
% xline(logf(20)) % check 20 Hz
% xline(logf(60)) % check line noise?
xticks([tick(2:7) logf(fStart+1) logf(fEnd+1)])
xticklabels([1:6 fStart fEnd])
xtickangle(0)

% Linear fit to log-log 1/f 
x = [log(meanF(fStart_idx)) log(meanF(fEnd_idx))];
slope = fitP(1); 
intercept = fitP(2); 
x = log(meanF(fStart_idx:fEnd_idx));
y = fitP(1)*x + fitP(2); 
rLine = refline(fitP(1),fitP(2));
rLine.Color = colors.green;
rLine.LineWidth = 1; 
rLine.LineStyle = '--';
rLine.DisplayName = '1/f extrapolation';
pl(3) = plot(x,y,'Color',colors.green,'LineWidth',1.5,'DisplayName','1/f fit'); 

if showN
    % --- Add n annotation ---
    nStr = sprintf('n = %d x 2 sessions',size(data.all.subject,2));
    nStrTxt = text(0.95*xl(2),yl(1)+0.1,nStr,'HorizontalAlignment','right','VerticalAlignment','bottom');
    nStrTxt.FontSize = 14;
    nStrTxt.FontName = 'Helvetica-Light';
end

if showLegend
    lgd = legend([pl(1:3) rLine]);
    lgd.FontSize = 10; 
    lgd.Box = 'off'; 
end

% get predicted log amps from log f
fOOF(1) = NaN; 
for i = 2:size(logf,2)
    fOOF(i) = slope*logf(i)+intercept;
end

pOOF(1) = NaN; 
for i = 2:size(f,2)
    pOOF(i) = invprctile( log10(dataBootGroup(i,:)) , fOOF(i) );
end

if saveFigs
    figTitle = sprintf('meg_manuscriptFigs_bootstrappedFFT_%s_%s',varOI,dateStr);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end
