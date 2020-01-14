function [tfAmps,freqoi,timeoi,fH] = meg_plotTF(data,selectedChannels,selectedFreq)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
%
% Karen Tian
% January 2020

%% args

if nargin<3
    selectedFreq = 8;
    disp('selectedChannels not specified, default 8')
end
if nargin<2 % default selects channels 1:3
    selectedChannels = 1:3;
    disp('selectedChannels not specified, default 1:3')
end
if nargin<1
    load('data.mat'); % load dummy data
end

%% setup
saveFigs = 0;

condNames = fieldnames(data);
nConds = numel(condNames);
nChannels = numel(selectedChannels);

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
eventTimes = [0 1000 1250 2100]; % check timing
eventNames = {'precue','T1','T2','response cue'};
tstart = -500; % -1000; 
tstop = 2800; % 2300;
t = tstart:tstop;

taper          = 'hanning';
foi            = 1:100;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfAmps = [];

% figures
ytick = 10:10:numel(foi);
xtick = 1:50:numel(toi);

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 
clims = [0 15];

nSamples = sz(1); 
nfft = 2^nextpow2(nSamples);
Fsample = 1000; 

%% time freq analysis

for iC = 1:nConds
    cond = condNames{iC}; 
    % erf.(thisField) = squeeze(nanmean(data.(thisField)(:,selectedChannels,:),2)); % time x trial
    
    y = fft(data.(cond)(:,selectedChannels,:),nfft)/nSamples; % freq x trial 
    f = Fsample/2*linspace(0,1,nfft/2+1); % frequencies
    amps = 2*abs(y(1:nfft/2+1,:)); % amp of freq x trial 
    ampsMean = nanmean(amps,3);
    
    % fieldtrip input
    vals = nanmean(data.(cond)(:,selectedChannels,:),3)'; % channels by samples
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(vals, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    
    specAmp = squeeze(mean(abs(spectrum),1)); % mean across channels
    tfAmps = specAmp'; % foi x time
    
    % store variables
    A.(cond).y = y; 
    A.(cond).f = f; 
    A.(cond).amps = amps; 
    A.(cond).ampsMean = ampsMean; 
    A.(cond).specAmp = specAmp;
    A.(cond).tfAmps = tfAmps;
    
end

%% spectrogram plot 

figure
set(gcf, 'Position',  [100, 100, 1000, 300*ceil(nConds/4)])
for iC = 1:nConds
    cond = condNames{iC};
    if nConds <= 4
        subplot (1,nConds,iC)
    end
    if nConds > 4
        subplot (ceil(nConds/3),nConds,iC)
    end
    hold on
    imagesc(A.(cond).tfAmps)
    xlim(xlims)
    ylim(ylims)
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    colorbar
    caxis(clims)
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    title(sprintf('%s \n channels %s ',und2space(cond), num2str(selectedChannels)))
end

%% normalized spectrogram plot 

figure
set(gcf, 'Position',  [100, 100, 1000, 300*ceil(nConds/4)])
for iC = 1:nConds
    cond = condNames{iC};
    if nConds <= 4
        subplot (1,nConds,iC)
    end
    if nConds > 4
        subplot (ceil(nConds/3),nConds,iC)
    end
    hold on
    for iF = 1:size(A.(cond).tfAmps,1)
        meanPow(iF) = nanmean(A.(cond).tfAmps(iF,:)); 
        for iT = 1:size(A.(cond).tfAmps,2)
            normPow(iF,iT) = A.(cond).tfAmps(iF,iT)/meanPow(iF)-1; 
        end
    end
    imagesc(normPow)
    xlim(xlims)
    ylim(ylims)
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    colorbar
    % caxis(clims)
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    title(sprintf('%s \n channels %s ',und2space(cond), num2str(selectedChannels)))
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

%% plot selected frequencies for all conditions

figure
for iH = 1:numel(selectedFreq)
    hold on
    freq = selectedFreq(iH);
    for iC = 1:nConds
        cond = condNames{iC};
        fVals = squeeze(A.(cond).tfAmps(freq,:));
        plot(toi*1000, fVals)
    end
    % xlim([t(1),t(end)])
    ylabel('amp')
    xlabel('time (ms)')
    vline(eventTimes,'k',eventNames)
    legend(condNames)
    title(sprintf('Frequency %d (Hz)',freq))
end

%% return figure handle

fH = sort(findobj('Type','figure'));

