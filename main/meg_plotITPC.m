function [A,fH,figNames] = meg_plotITPC(data,sessionDir,p,selectedChannels,selectedFreq)

% MEG_PLOTITPC(data,sessionDir,p,selectedChannels,selectedFreq)
%
% INPUT
% data
%   structure of condition cells containing
%   data matrix, time x channels x trials
% sessionDir 
%   eg: 'R0817_20171212'
% selectedChannels
%   vector of meg channels of interest 
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
% July 2020

%% args

if nargin<4
    selectedFreq = 20;
    disp('selectedFreqs not specified, default 20 Hz')
end
if nargin<3 % default selects channels 1:3
    selectedChannels = [20,23,36,43,60];
    disp('selectedChannels not specified, default [20,23,36,43,60]')
end
if nargin<2
    p = meg_params('TANoise_ITPC'); 
    disp('timing parameters not specified, using TANoise ITPC')
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
figNames = []; 
plotFigs = 1; % 0  

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
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% params 
% timing 
t = p.tstart:p.tstop;

taper          = 'hanning';
foi            = selectedFreq; % 100;
t_ftimwin      = 10 ./ foi;
toi            = p.tstart/1000:0.01:p.tstop/1000;
tfAmps = []; tfAmpsAvg = [];
tfPows = []; tfPowsAvg = [];

padTotal = ceil(p.trialTime/p.fSample); 
padPre = ceil((padTotal*p.fSample-p.trialTime)/2); 
padPost = floor((padTotal*p.fSample-p.trialTime)/2); 
toiPad = (p.tstart-padPre)/1000:0.01:(p.tstop+padPost)/1000; 
tPad = p.tstart-padPre:p.tstop+padPost; 
xtick = 1:80:numel(toiPad);

% figures
ytick = 10:10:numel(foi);
% xtick = 1:50:numel(toi);
% xtickms = p.tstart:500:p.tstop;
xlims = [size(toi,1),size(toi,2)]; 

% xlimsPad = [1 padTotal*100]; 

nSamples = sz(1); 
% nfft = 1000; % 2^nextpow2(nSamples); % 1000; 
Fsample = p.fSample; 
width = 8; 

%% time freq analysis
%     y = fft(data.(cond)(1:p.prestim*1000,selectedChannels,:),nfft)/nSamples; % freq x trial 
%     f = Fsample/2*linspace(0,1,nfft/2+1); % frequencies
%     amps = 2*abs(y(1:nfft/2+1,:,:)); % amp of freq x trial 
%     ampsMean = nanmean(amps,3);
    
%% fieldtrip wavelet convolution 
for iC = 1:numel(condNames)
    for iCh = 1:nChannels
        % single trial
        meanSpectrumAll = [];
        vals = [];
        vals = squeeze(nanmean(data.(condNames{iC})(:,selectedChannels(iCh),:),2))'; % trials by samples
        % vals = vals(:,1:7000); % November 3 ITPC spectrogram edit
        [spectrum,freqoi,timeoi] = ft_specest_wavelet(vals, t/1000, ...
            'freqoi', foi, 'width', width);

        tfAmpsTrials = abs(spectrum); 
        tfPowsTrials = tfAmpsTrials.^2; 

        meanSpectrum = squeeze(mean(abs(spectrum),1,'omitnan')); % t x f for one channel, average trials
        meanSpectrumAll = cat(3,meanSpectrumAll,meanSpectrum); % t x f x channel spectrum
        itpc = squeeze(abs(mean(exp(1i*angle(spectrum)),1,'omitnan')));
        tfAmps = squeeze(permute(meanSpectrumAll,[3,2,1])); % foi x time
        tfPows = tfAmps.^2;
        
        % trial average
        %         meanSpectrumAllAvg = [];
        %         valsAvg = [];
        %         valsAvg = nanmean(vals,1);
        %         [spectrumAvg,freqoiAvg,timeoiAvg] = ft_specest_wavelet(valsAvg, t/1000, ...
        %             'freqoi', foi, 'width', width);
        %         meanSpectrumAvg = squeeze(nanmean(abs(spectrumAvg),1)); % t x f for one channel
        %         meanSpectrumAllAvg = cat(3,meanSpectrumAllAvg,meanSpectrumAvg); % t x f x channel spectrum
        %         tfAmpsAvg = squeeze(permute(meanSpectrumAllAvg,[3,2,1])); % foi x time
        %         tfPowsAvg = tfAmpsAvg.^2;
        
        A.(condNames{iC}).selectedChannels = selectedChannels;
        %         A.(condNames{iC}).tfAmps(:,:,iCh) = tfAmps; % f x t x  ch
        %         A.(condNames{iC}).tfPows(:,:,iCh) = tfPows; % f x t x  ch
        A.(condNames{iC}).tfAmpsTrials(:,:,:,iCh) = tfAmpsTrials; 
        A.(condNames{iC}).tfPowsTrials(:,:,:,iCh) = tfPowsTrials; 
        A.(condNames{iC}).tfAmps(:,:) = mean(tfAmps,3,'omitnan'); % f x t x  ch
        A.(condNames{iC}).tfPows(:,:) = mean(tfPows,3,'omitnan'); % f x t x  ch
        A.(condNames{iC}).ITPC(:,:,iCh) = itpc; % f x t x  ch
        A.(condNames{iC}).spectrum(:,:,iCh) = squeeze(spectrum); % complex values  
        A.(condNames{iC}).phaseAngle(:,iCh,:,:) = angle(spectrum);
    end
    A.(condNames{iC}).ITPCMean = nanmean(A.(condNames{iC}).ITPC,3); % average channels 
    % A.(condNames{iC}).phaseAngleCh = squeeze(nanmean(A.(condNames{iC}).phaseAngle,1));
end

if plotFigs 
    %% plot time series check 
    figure
    set(gcf, 'Position',  [100, 100, 800, 300])
    hold on 
    for iC = 1:nConds
        vals = []; 
        vals = nanmean(data.(condNames{iC})(:,selectedChannels,:),3); 
        vals = nanmean(vals,2); % average channels 
        plot(1:size(vals,1),vals,'Color',p.cueColors(iC,:),'LineWidth',2)
    end
    vline(p.eventTimes,':k')
    set(gca,'TickDir','out');
    ax = gca;
    ax.LineWidth = 1.5;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize = 12;
    xlim([p.tstart,p.tstop])
    ylabel('Amplitude (T)')
    xlabel('Time (ms)')
    title(sprintf('%s Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))
    legend(condNames) 

    %% plot 20Hz ITPC 
    figure
    set(gcf, 'Position',  [100, 100, 800, 300])
    hold on 
    for iC = 1:nConds
        vals = []; 
        vals = A.(condNames{iC}).ITPCMean; % average channels 
        plot(t,vals,'Color',p.cueColors(iC,:),'LineWidth',2)
    end
    vline(p.eventTimes,':k')
    set(gca,'TickDir','out');
    ax = gca;
    ax.LineWidth = 1.5;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize = 12;
    xlim([p.tstart,p.tstop])
    ylabel('ITPC')
    xlabel('Time (ms)')
    title(sprintf('%s Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))
    legend(condNames) 

    %% plot single trial power 

    % figure
    % set(gcf, 'Position',  [100, 100, 800, 300])
    % hold on 
    % for iC = 1:nConds
    %     errBar = shadedErrorBar(t,nanmean(A.(condNames{iC}).tfPows,2),nanstd(A.(condNames{iC}).tfPows,[],2)/sqrt(nChannels));
    %     errBar.patch.FaceColor = p.cueColors(iC,:); 
    %     errBar.edge(1).Color = p.cueColors(iC,:);
    %     errBar.edge(2).Color = p.cueColors(iC,:);
    %     plot(t,nanmean(A.(condNames{iC}).tfPows,2),'Color',p.cueColors(iC,:),'LineWidth',2)
    % end
    % vline(p.eventTimes,':k')
    % set(gca,'TickDir','out');
    % ax = gca;
    % ax.LineWidth = 1.5;
    % ax.XColor = 'black';
    % ax.YColor = 'black';
    % ax.FontSize = 12;
    % ylabel('Single Trial 20 Hz Power (T^2)')
    % xlabel('Time (ms)')
    % title(sprintf('%s Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))

    %% plot average trial power 

    % figure
    % set(gcf, 'Position',  [100, 100, 800, 300])
    % hold on 
    % for iC = 1:nConds
    %     plot(t,nanmean(A.(condNames{iC}).tfPowsAvg,2),'Color',p.cueColors(iC,:),'LineWidth',2)
    % end
    % vline(p.eventTimes,':k')
    % set(gca,'TickDir','out');
    % ax = gca;
    % ax.LineWidth = 1.5;
    % ax.XColor = 'black';
    % ax.YColor = 'black';
    % ax.FontSize = 12;
    % ylabel('Average Trial 20 Hz Power (T^2)')
    % xlabel('Time (ms)')
    % title(sprintf('%s Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))

    %% plot ITPC 

    % figure
    % set(gcf, 'Position',  [100, 100, 800, 300])
    % hold on 
    % for iC = 1:nConds
    %     errBar = shadedErrorBar(t,A.(condNames{iC}).ITPCMean,nanstd(A.(condNames{iC}).ITPC,[],2)/sqrt(nChannels));
    %     errBar.patch.FaceColor = p.cueColors(iC,:); 
    %     errBar.edge(1).Color = p.cueColors(iC,:);
    %     errBar.edge(2).Color = p.cueColors(iC,:);
    %     plot(t,A.(condNames{iC}).ITPCMean,'Color',p.cueColors(iC,:),'LineWidth',2)
    % end
    % vline(p.eventTimes,':k')
    % set(gca,'TickDir','out');
    % ax = gca;
    % ax.LineWidth = 1.5;
    % ax.XColor = 'black';
    % ax.YColor = 'black';
    % ax.FontSize = 12;
    % ylabel('ITPC')
    % xlabel('Time (ms)')
    % title(sprintf('%s Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))

    %% ITPC phase angle at T1 (average trials per channel)

    % figure
    % set(gcf, 'Position',  [100, 100, 300, 400])
    % for iC = 1:nConds
    %     theta = A.(condNames{iC}).phaseAngleCh(:,p.eventTimes(2)+abs(p.tstart))';
    %     for i = 1:nChannels
    %         pp = polarplot([theta(i),theta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',0.5);
    %         pp.Color(4) = 0.5; 
    %         hold on
    %     end
    %     pp = polarplot([circ_mean(theta,[],2),nanmean(theta,2)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',5);
    %     pp.Color(4) = 0.7; 
    % end
    % grid off 
    % ppAx = gca;
    % ppAx.RColor = [0 0 0];
    % ppAx.LineWidth = 1; 
    % % ppAx.ThetaTick = 0:30:330; 
    % % ppAx.ThetaTick = []; 
    % ppAx.RTick = [0,1];
    % rlim([0,1])
    % title(sprintf('Phase Angle at T1 \n %s \n Channels %s',und2space(sessionDir),num2str(selectedChannels(1:nChannels))))

    %% ITPC phase angle at T1 (all trials, top channel)

    % figure
    % set(gcf, 'Position',  [100, 100, 300, 400])
    % for iC = 1:nConds
    %     theta = A.(condNames{iC}).phaseAngle(:,1,p.eventTimes(2)+abs(p.tstart))';
    %     theta = theta(~isnan(theta)); 
    %     for i = 1:numel(theta)
    %         pp = polarplot([theta(i),theta(i)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',0.5);
    %         pp.Color(4) = 0.5; 
    %         hold on
    %     end
    %     pp = polarplot([circ_mean(theta,[],2),nanmean(theta,2)],[0,1],'Color',p.cueColors(iC,:),'LineWidth',6);
    %     pp.Color(4) = 0.6; 
    % end
    % grid off 
    % ppAx = gca;
    % ppAx.RColor = [0 0 0];
    % ppAx.LineWidth = 1; 
    % % ppAx.ThetaTick = 0:30:330; 
    % % ppAx.ThetaTick = []; 
    % ppAx.RTick = [0,1];
    % rlim([0,1])
    % title(sprintf('Phase Angle at T1 \n %s \n Channels %s',und2space(sessionDir),num2str(selectedChannels(1:1))))
end
%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
% figNames = {'timeSeries','singleTrialPower','averageTrialPower','ITPC','phaseAngleT1','phaseAngleT1Ch1'};
% figNames = {'timeSeries','singleTrialPower','ITPC','phaseAngleT1'};
figNames = {'timeSeries','ITPC20Hz'};

