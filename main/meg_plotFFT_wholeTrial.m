function [A,fH,figNames] = meg_plotFFT_wholeTrial(data,p,selectedChannels,selectedFreq)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% INPUT
% data
%   structure of condition cells containing
%   data matrix, time x channels x trials
% selectedChannels
%   vector of meg channels to inspect 
% selectedFreq  
%   vector of frequencies of interest
% p 
%   parameters (for expt short type)
% OUPUT
% A
%   structure of time freq amp
% fH
%   figure handle
% figNames 
%   figure names 
%
% Karen Tian
% January 2020

%% args

if nargin<4
    selectedFreq = 8:14;
    disp('selectedFreqs not specified, default 8-14 Hz (alpha band (Foxe 2011))')
end
if nargin<3 % default selects channels 1:3
    selectedChannels = [20,23,36,43,60];
    disp('selectedChannels not specified, default [20,23,36,43,60]')
end
if nargin<2
    p = meg_params('TA2_Analysis'); 
    disp('timing parameters not specified, using TA2 analysis')
end
if nargin<1
    load('D3.mat'); % load dummy data
    data = D3; 
end

%% setup
condNames = fieldnames(data);
nConds = numel(condNames);
nChannels = numel(selectedChannels);
nFreqs = numel(selectedFreq); 

% index selected channels; 
data.(condNames{1}) = data.(condNames{1})(:,selectedChannels,:);  

% chop off last time pt? for whole trial 
% data.(condNames{1})(end,:,:) = [];

% chop off all post precue
% data.(condNames{1})(end-5000:end,:,:) = [];

% pretarget only 
data.(condNames{1})(end-3950:end,:,:) = [];

%% checks

% check data
sz = size(data.(condNames{1}));
if numel(sz)<2 || numel(sz)>3
    error('data is expected to be at least 2-dimensional, with time as the first dimension, and channels second')
elseif numel(sz)==2 
    nTrials = 1; 
elseif numel(sz)==3
    nTrials = sz(3); 
end
nAllChannels = sz(2); 
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% params 

% timing 
t = p.tstart:p.tstop;

taper          = 'hanning';
foi            = 1:100;
t_ftimwin      = 10 ./ foi;
toi            = p.tstart/1000:0.01:p.tstop/1000;
tfAmps = [];

% figures
ytick = 10:10:numel(foi);
xtick = 1:50:numel(toi);
xtickms = p.tstart:500:p.tstop;

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 

nSamples = sz(1); 
nfft = 1000; % nSamples/10; % try for integer freqs? 2^nextpow2(nSamples);
Fsample = 1000; 

cueColors = [122/255 142/255 194/255; 225/255 124/255 96/255; 128/255 128/255 128/255];  % cueT1, cueT2, neutral 
colorAlpha = 0.75; % transparency for plots 

%% time freq analysis

for iC = 1:nConds
    cond = condNames{iC}; 
    
    y = fft(data.(cond),nfft)/nSamples; % time x freq x trial 
    f = Fsample/2*linspace(0,1,nfft/2+1); % frequencies
    amps = 2*abs(y(1:nfft/2+1,:,:)); % amp of freq x trial 
    ampsMean = nanmean(amps,3);
    
    % fieldtrip wavelet convolution
    meanSpectrumAll = []; 
    for iCh = 1:numel(selectedChannels)
        vals = squeeze(data.(cond)(:,iCh,:))'; % trials by samples
        % valspad = padarray(vals,[0 5*p.fSample-nSamples],'post'); % pad
        % so integer time freqs
        [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(vals, t/1000, ...
            'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
            'taper', taper, 'dimord', 'chan_time_freqtap');
        meanSpectrum = squeeze(nanmean(abs(spectrum),1)); % t x f for one channel 
        meanSpectrumAll = cat(3,meanSpectrumAll,meanSpectrum); % t x f x channel spectrum 
    end
    
    tfAmps = permute(meanSpectrumAll,[3,2,1]); % foi x time
    meanTfAmps = squeeze(nanmean(tfAmps,1)); 
    tfPows = tfAmps.^2; 
    meanTfPows = meanTfAmps.^2; 
    
    % normalize
    meanMeanAmps = nanmean(meanTfAmps,2);
    meanMeanPows = nanmean(meanTfPows,2);
    normAmps = meanTfAmps./meanMeanAmps-1;
    normPows = meanTfPows./meanMeanPows-1;
    
    % store variables
    A.(cond).y = y; 
    A.(cond).f = f; 
    A.(cond).amps = amps; 
    A.(cond).ampsMean = ampsMean; 
    A.(cond).tfAmps = tfAmps;
    A.(cond).meanTfAmps = meanTfAmps;
    A.(cond).tfPows = tfPows; 
    A.(cond).meanTfPows = meanTfPows;
    A.(cond).normAmps = normAmps; 
    A.(cond).normPows = normPows; 
    A.(cond).selectedChannels = selectedChannels; 
end

%% spectrogram plot 

figure
for iC = 1:nConds
    cond = condNames{iC};
    if nConds <= 4
        subplot (1,nConds,iC)
        set(gcf, 'Position',  [100, 100, 300*nConds, 300*ceil(nConds/4)])
    end
    if nConds > 4
        subplot (ceil(nConds/3),nConds,iC)
        set(gcf, 'Position',  [100, 100, 300*ceil(nConds/4), 300*ceil(nConds/4)])
    end
    hold on
    imagesc(A.(cond).meanTfPows) 
    xlim(xlims)
    ylim(ylims)
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    colorbar
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    title(sprintf('%s \n channels %s ',und2space(cond), num2str(selectedChannels)))
end

%% normalized spectrogram plot 

figure
set(gcf, 'Position',  [100, 100, 300*ceil(nConds/4), 300*ceil(nConds/4)])
for iC = 1:nConds
    cond = condNames{iC};
    if nConds <= 4
        set(gcf, 'Position',  [100, 100, 300*nConds, 300*ceil(nConds/4)])
        subplot (1,nConds,iC)
    end
    if nConds > 4
        set(gcf, 'Position',  [100, 100, 300*ceil(nConds/4), 300*ceil(nConds/4)])
        subplot (ceil(nConds/3),nConds,iC)
    end
    hold on
    imagesc(A.(cond).normPows)
    xlim(xlims)
    ylim(ylims)
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    colorbar
    % caxis([min(min(A.(cond).normPows)) max(max(A.(cond).normPows))])
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    title(sprintf('power normalized %s \n channels %s ',und2space(cond), num2str(selectedChannels)))
end

%% fourier plot
figure
for iC = 1:nConds
    cond = condNames{iC};
    % set(gcf, 'Position',  [100, 100, 800, 300])
    subplot (nConds,1,iC)
    loglog(A.(cond).f, A.(cond).ampsMean) 
    hold on
    selectChAmps = nanmean(A.(cond).ampsMean(:,:),2); 
    loglog(A.(cond).f, selectChAmps,'k','LineWidth',3) 
    xlim([A.(cond).f(1) A.(cond).f(end)])
    % title(sprintf('cond %s, whole trial',und2space(cond)))
    title(sprintf('cond %s, pre target',und2space(cond)))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    
    set(gca,'TickDir','out');
    ax = gca;
    box off
    ax.LineWidth = 1.5;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize = 12;
end

%% fourier plot average

figure
set(gcf, 'Position',  [100, 100, 800, 300])
for iC = 1:nConds
    cond = condNames{iC};
    selectChAmps = nanmean(A.(cond).ampsMean(:,:),2); 
    figFourier = loglog(A.(cond).f, selectChAmps,'Color',cueColors(iC,:));
    % figFourier.Color(4) = colorAlpha; 
    hold on
    xlim([A.(cond).f(1) A.(cond).f(end)])
    title('fft')
    legend(condNames)
    xlabel('frequency (Hz)')
    ylabel('amp')
end

%% plot selected frequency range for all conditions

figure
hold on
for iC = 1:nConds
    cond = condNames{iC};
    fVals = squeeze(nanmean(A.(cond).meanTfAmps(selectedFreq,:),1));
    figfreq = plot(toi, fVals,'Color',cueColors(iC,:),'LineWidth',2);
    figfreq.Color(4) = colorAlpha;   
end
xlim([toi(1),toi(end)])
ylabel('amp')
xlabel('time (s)')
vline(p.eventTimes/1000,'k',p.eventNames)
legend(condNames)
title(sprintf('Frequency %d - %d (Hz)',selectedFreq(1),selectedFreq(end)))

%% diff TF 
% condition1 = A.cueT1.normPows; 
% condition2 = A.cueT2.normPows; 
% diffc2c1 = condition2 - condition1; % difference
% % meanDiff = nanmean(diffc2c1,2); % mean diff pow 
% % diff = diffc2c1./meanDiff-1; % then normalize 
% 
% figure 
% set(gcf,'Position',[100 100 300 350])
% hold on
% imagesc(diffc2c1)
% xlim(xlims)
% ylim(ylims)
% xlabel('time (s)')
% ylabel('frequency (Hz)')
% colorbar
% meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
% title('cue T2 - cue T1 (power, norm)')

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
% figNames = {'TF','TFnorm','FFTAllCh','FFTSelectCh','Freq','TF_diffT2T1'}; 
figNames = {'TF','TFnorm','FFTAllCh','FFTSelectCh','Freq'}; 

