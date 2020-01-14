function fH = meg_plotTF(data,selectedChannels)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
%
% Karen Tian
% January 2020

%% args

if nargin<2 % default selects channels 1:3
    selectedChannels = 1:3;
    disp('selectedChannels not specified, default 1:3')
end
if nargin<1
    load('data.mat'); % load dummy data
end

%% setup
saveFigs = 0;

fieldName = fieldnames(data);
nFields = numel(fieldName);
nChannels = numel(selectedChannels);

% timing 
eventTimes = [0 1000 1250 2100]; % check timing
eventNames = {'precue','T1','T2','response cue'};
tstart = -1000; % for targets
tstop = 2300;
t = tstart:tstop;

%% checks

% check data
sz = size(data.(fieldName{1}));
if numel(sz)~=3
    error('data is expected to be 3-dimensional, with trials as the last dimension')
end
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% ft time freq params 

taper          = 'hanning';
foi            = 1:100;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfAmps = [];

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 
clims = [0 15];

nSamples = sz(1); 
nfft = 2^nextpow2(nSamples);

%% FFT on high SNR channels

for iF = 1:nFields
    thisField = fieldName{iF}; 
    targetERF.(thisField) = squeeze(nanmean(data.(thisField)(:,selectedChannels,:),2)); % time x trial 
    targetY.(thisField) = fft(targetERF.(thisField),nfft)/nSamples; % freq x trial 
    targetAmps.(thisField) = 2*abs(targetY.(thisField)(1:nfft/2+1,:)); % amp of freq x trial 
    
    vals.(thisField) = nanmean(data.(thisField)(:,selectedChannels,:),3)'; % channels by samples
    
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(vals.(thisField), t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    
    specAmp = squeeze(mean(abs(spectrum),1)); % mean across channels
    tfAmps.(thisField) = specAmp';
   
    figure
    hold on 
    imagesc(tfAmps.(thisField))
    xlim([1,301])
    ylim([1,101])
    meg_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes(1:3));
    title(sprintf('%s',und2space(thisField)))
end
 
%% OTHER JUNK TESTING
%%%%%%%%%%%%%%%%%%%%%

% ssvefFreq = 30;
% noiseDist = 5;
% freqIdx = find(foi==ssvefFreq);
% signal = tfAmps(freqIdx,:,:);
% noise = tfAmps(freqIdx + [-noiseDist noiseDist],:,:);
% 
% % figures
% ytick = 10:10:numel(foi);
% xtick = 51:50:numel(toi);
% clims = [0 15];
% diffClims = [-10 10];
% eventTimes = 0
% 
% fH(2) = figure;
% for iTrig = 1:nTrigs
%     subplot(1,2,iTrig)
%     imagesc(tfAmps(:,:,iTrig),clims)
%     rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
%     if iTrig==nTrigs
%         xlabel('time (ms)')
%         ylabel('frequency (Hz)')
%     end
%     title(trigNames{iTrig})
% end
% rd_supertitle(['channel' sprintf(' %d', channels)]);
% rd_raiseAxis(gca);
% 
% 
% 
% 
% 
% %% ft read data
% 
% sqdfile = '/Users/kantian/Dropbox/Data/TA2/MEG/R0817_20181120/preproc/R0817_TA2_11.20.18_run01_bietfp.sqd';
% 
% dat = ft_read_data(sqdfile);
% hdr = ft_read_header(sqdfile);
% 
% cfg                     = [];
% cfg.dataset             = sqdfile;
% cfg.trialdef.prestim    = 0.2; %0; %0.5; % sec
% cfg.trialdef.poststim   = 2.3; %2.3; %3.1;
% cfg.trialdef.trig       = [168,167]; %[168,167]; %[161:164,167]; %168 = precue, 167=blank
% threshold = 2.5;
% [trl,Events]            = mytrialfun_all(cfg,threshold,[]);
% 
% data          = ft_preprocessing(struct('dataset',sqdfile,...
%     'channel','MEG','continuous','yes','trl',trl));
% 
% %% ft time freq analysis 
% 
% cfg              = [];
% cfg.output       = 'pow';
% cfg.channel      = 'MEG';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
% cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% TFRhann = ft_freqanalysis(cfg, data);
% 
% %% plot
% 
% cfg              = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'absolute';
% cfg.maskstyle    = 'saturation';
% cfg.zlim         = [-3e-27 3e-27];
% cfg.channel      = 'MRC15';
% cfg.interactive  = 'no';
% figure
% ft_singleplotTFR(cfg, TFRhann);
% 
% %%
% 
% power = TFRhann.powspctrm; % channels x freqs x trials 
% powerCh = power(selectedChannels, :, :); 
% meanPowerCh = nanmean(powerCh,3); 
% 
% imagesc(meanPowerCh)
% 
% 
% %% return figure handle
% 
% rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes)

