function [A,fH,figNames] = meg_plotTF(data,selectedChannels,selectedFreq)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% INPUT
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
% selectedChannels
%   vector of meg channels to inspect 
% selectedFreq  
%   vector of frequencies of interest
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

if nargin<3
    selectedFreq = 8:14;
    disp('selectedFreqs not specified, default 8-14 Hz (alpha band (Foxe 2011))')
end
if nargin<2 % default selects channels 1:3
    selectedChannels = [20,23,36,43,60];
    disp('selectedChannels not specified, default [20,23,36,43,60]')
end
if nargin<1
    load('data.mat'); % load dummy data
end

%% setup
saveFigs = 0;

condNames = fieldnames(data);
nConds = numel(condNames);
nChannels = numel(selectedChannels);
nFreqs = numel(selectedFreq); 

%% checks

% check data
sz = size(data.(condNames{1}));
if numel(sz)~=3
    error('data is expected to be 3-dimensional, with trials as the last dimension')
end
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% params 

% timing 
p = meg_params('TA2');
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
clims = [0 15];

nSamples = sz(1); 
nfft = 2^nextpow2(nSamples);
Fsample = 1000; 

%% time freq analysis

for iC = 1:nConds
    cond = condNames{iC}; 
    
    y = fft(data.(cond),nfft)/nSamples; % freq x trial 
    f = Fsample/2*linspace(0,1,nfft/2+1); % frequencies
    amps = 2*abs(y(1:nfft/2+1,:,:)); % amp of freq x trial 
    ampsMean = nanmean(amps,3);
    
    % fieldtrip wavelet convolution
    meanSpectrumAll = []; 
    for iCh = selectedChannels
        vals = squeeze(data.(cond)(:,iCh,:))'; % trials by samples
        [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(vals, t/1000, ...
            'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
            'taper', taper, 'dimord', 'chan_time_freqtap');
        meanSpectrum = squeeze(nanmean(abs(spectrum),1)); % t x f for one channel 
        meanSpectrumAll = cat(3,meanSpectrumAll,meanSpectrum); % t x f x channel spectrum 
    end
    
    tfAmps = permute(meanSpectrumAll,[3,2,1]); % foi x time
    meanTfAmps = squeeze(nanmean(tfAmps,1)); 
    
    % store variables
    A.(cond).y = y; 
    A.(cond).f = f; 
    A.(cond).amps = amps; 
    A.(cond).ampsMean = ampsMean; 
    A.(cond).tfAmps = tfAmps;
    A.(cond).meanTfAmps = meanTfAmps;
    
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
    imagesc(A.(cond).meanTfAmps) 
    xlim(xlims)
    ylim(ylims)
    % xticklabels(num2str(xtickms))
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
    for iF = 1:numel(foi)
        meanPow(iF) = nanmean(A.(cond).meanTfAmps(iF,:)); 
        for iT = 1:numel(toi)
            normPow(iF,iT) = A.(cond).meanTfAmps(iF,iT)/meanPow(iF)-1; 
        end
    end
    imagesc(normPow)
    xlim(xlims)
    ylim(ylims)
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    colorbar
    caxis([-0.2 0.2])
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    title(sprintf('normalized %s \n channels %s ',und2space(cond), num2str(selectedChannels)))
end

%% fourier plot

figure
for iC = 1:nConds
    cond = condNames{iC};
    % set(gcf, 'Position',  [100, 100, 800, 300])
    subplot (nConds,1,iC)
    loglog(A.(cond).f, A.(cond).ampsMean) 
    hold on
    selectChAmps = nanmean(A.(cond).ampsMean(:,selectedChannels),2); 
    loglog(A.(cond).f, selectChAmps,'k','LineWidth',3) 
    xlim([A.(cond).f(1) A.(cond).f(end)])
    title(sprintf('%s',und2space(cond)))
    xlabel('frequency (Hz)')
    ylabel('amp')
end

%% fourier plot average

figure
set(gcf, 'Position',  [100, 100, 800, 300])
for iC = 1:nConds
    cond = condNames{iC};
    selectChAmps = nanmean(A.(cond).ampsMean(:,selectedChannels),2); 
    loglog(A.(cond).f, selectChAmps) 
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
    plot(toi, fVals)
end
xlim([toi(1),toi(end)])
ylabel('amp')
xlabel('time (ms)')
vline(p.eventTimes/1000,'k',p.eventNames)
legend(condNames)
title(sprintf('Frequency %d - %d (Hz)',selectedFreq(1),selectedFreq(end)))

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
figNames = {'TF','TFnorm','FFTAllCh','FFTSelectCh','Freq'}; 

