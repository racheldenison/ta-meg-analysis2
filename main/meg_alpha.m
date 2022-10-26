function [A,fH,figNames] = meg_alpha(eyesClosedFile,sessionDir)

% MEG_ALPHA(eyesClosedFile,sessionDir)
%
% INPUTS
% eyesClosedFile
%   sqd
% sessionDir
%   eg: 'R1507_20190627'
% 
% OUPUT
% fH
%   figure handle 
%   topo plot of individual peak alpha power ± 2Hz
% figNames 
% alpha 
%   peak individual alpha (Hz)
%
% Karen Tian
% January 2020

%% args
if nargin<1 
    disp('enter eyes closed .sqd file path')
end

%% read data 
data = sqdread(eyesClosedFile); % read data 
p = meg_params('TANoise_Analysis'); % get parameters 
data = data(:,p.megChannels); % read only meg data

%% setup
idx = 5000:14999; % take middle chunk of data 
data = data(idx,:); 
nSamples = size(data,1); 
nfft = 1000; 

% timing 
sz = size(data);
t = 1:sz(1); 

taper          = 'hanning';
foi            = 1:25;
t_ftimwin      = 10 ./ foi;
toi            = t(1)/p.fSample:0.01:t(end)/p.fSample;
tfAmps = [];
window = 8:14; % Hz, to look for peak alpha, Engel 2019

% figures
ytick = 10:10:numel(foi);
xtick = 1:50:numel(toi);

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 

load('data_hdr.mat')

% nSamples = t(end); 
% nfft = 2^nextpow2(t(end));
% nfft = 50; 

%% fft
y = fft(data,nfft)/nSamples; % freq x trial
f = p.fSample/2*linspace(0,1,nfft/2+1); % frequencies
amps = 2*abs(y(1:nfft/2+1,:,:)); % amp of freq x trial
ampsMean = nanmean(amps,3);

%% fieldtrip wavelet convolution time freq analysis 
meanSpectrumAll = [];
for iCh = p.megChannels
    vals = squeeze(data(:,iCh))'; % channels by samples
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(vals, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    meanSpectrum = squeeze(nanmean(abs(spectrum),1)); % t x f for one channel
    meanSpectrumAll = cat(3,meanSpectrumAll,meanSpectrum); % t x f x channel spectrum
end

tfAmps = permute(meanSpectrumAll,[3,2,1]); % foi x time
meanTfAmps = squeeze(nanmean(tfAmps,1)); % average across channels 

% store variables
A.y = y;
A.f = f;
A.amps = amps;
A.ampsMean = ampsMean;
A.tfAmps = tfAmps;
A.meanTfAmps = meanTfAmps;
A.tfPows = tfAmps.^2;
A.meanTfPows = meanTfAmps.^2;

%% spectrogram plot 
figure
hold on
imagesc(A.meanTfPows)
xlim(xlims)
ylim(ylims)
xlabel('Time (samples) ')
ylabel('Frequency (Hz)')
colorbar
h = colorbar;
set(get(h,'title'),'string','Power');
title(sprintf('Eyes Closed %s ', und2space(sessionDir)))

%% tf output alpha peak, avg all channels -> find peak 
meanmeanTfAmps = nanmean(A.meanTfAmps,2); 
[peak,peakIdx] = max(meanmeanTfAmps(window)); 
A.tfAlpha = window(peakIdx); 

figure
plot(meanmeanTfAmps)
hold on
xlim([window(1) window(end)])
vline(A.tfAlpha,'k')
title(sprintf('tf alpha %d Hz %s',A.tfAlpha,und2space(sessionDir)))

%% tf alpha peak, single channels -> average peak 
A.chTfPows = nanmean(A.tfAmps,3); % ch x foi; average across time
for iCh = 1:size(A.chTfPows,1) % find max tf pow within within 
    [peak(iCh),peakIdx(iCh)] = max(A.chTfPows(iCh,window)); 
end

A.tfChAlpha = window(peakIdx); 

figure
plot(A.chTfPows)
hold on
xlim([window(1) window(end)])
vline(A.tfChAlpha,'k')
title(sprintf('Single channel tf alpha %d Hz %s',A.tfAlpha,und2space(sessionDir)))

%% fourier, avg all channels 
ampsMeanMean = nanmean(A.ampsMean,2);

windowIdx = find(A.f>window(1) & A.f<window(end));
[peak,peakIdx] = max(ampsMeanMean(windowIdx)); 

figure
set(gcf, 'Position',  [100, 100, 800, 300])
loglog(A.f, ampsMeanMean)
hold on
xlim([A.f(1) A.f(end)])
title('fft')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% plot peak frequency within window
A.alpha = A.f(peakIdx)+window(1); 
vline(A.alpha,'k')
title(sprintf('%s alpha %.2f Hz',und2space(sessionDir),A.alpha))

%% fourier, single channels 
windowIdx = find(A.f>window(1) & A.f<window(end)); % find frequency indices within alpha window 
for iCh = 1:size(A.ampsMean,2) % for each channel, find max and max index 
    [peak(iCh),peakIdx(iCh)] = max(A.ampsMean(windowIdx,iCh));
end

% do a mov mean? to smooth fft and reduce frequencies?
% f10 = movmean(A.f,10); 
% windowIdx = find(f10>window(1) & f10<window(end));
% 
% for iCh = 1:size(A.ampsMean,2) % for each channel, find max and max index 
%     [peak(iCh),peakIdx(iCh)] = max(A.ampsMean(windowIdx,iCh));
% end

for iCh = 1:size(A.ampsMean,2) % for each channel, find max and max index 
    [peak(iCh),peakIdx(iCh)] = max(A.ampsMean(windowIdx,iCh));
end

A.fftChAlpha = A.f(peakIdx)+window(1); % store 

% fDownsample = movmean(A.f,1000); 
% ampsDownsample = movmean(A.ampsMean,1000); 

figure
set(gcf, 'Position',  [100, 100, 800, 300])
loglog(A.f', A.ampsMean)
% loglog(fDownsample, ampsDownsample)
hold on
xlim([A.f(1) A.f(end)])
title('fft')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% plot peak frequency within window
vline(A.fftChAlpha,'k')
title(sprintf('Single channel fft peaks %s alpha %.2f Hz',und2space(sessionDir),A.alpha))

%% histogram of fft single ch peaks 
edges = window(1):0.25:window(end); 
[Npeaks,edgesPeak] = histcounts(A.fftChAlpha,edges); % count peaks within bin size 

% valsPeaks = []; 
% for i = 1:numel(edges)-1
%     replicateVals = Npeaks(i); 
%     vals = repmat(replicateVals,1,replicateVals); 
%     valsPeaks = [valsPeaks vals]; 
% end
% yVals = ones(1,size(valsPeaks,2)); 
% figure % beeswarm plot 
% beeswarm(valsPeaks',yVals,'sort_style','up','corral_style','random','dot_size',3,'overlay_style','ci');

A.fftChAlphaAvg = nanmean(A.fftChAlpha); 
figure
histogram('BinEdges',edgesPeak,'BinCounts',Npeaks)
title(sprintf('%s single ch fft peaks distribution, avg = %.2f Hz',und2space(sessionDir),A.fftChAlphaAvg))
ylabel('Count (number of channels)')
xlabel('Frequency (0.25Hz bins)')
set(gca,'TickDir','out');
vline(A.fftChAlphaAvg,':k')
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off
set(gca,'TitleFontSizeMultiplier',1.25)

%% topo of channel peak alpha
figure
meg_topoplot(A.fftChAlpha')
caxis([window(1) window(end)])
colorbar
title(sprintf('%s peak alpha power per channel, 8-14 Hz',und2space(sessionDir)))

%% field trip multiploter tester
% cfg                     = [];
% cfg.dataset             = eyesClosedFile;
% eyesClosedData          = ft_preprocessing(struct('dataset',eyesClosedFile,...
%     'channel','MEG','continuous','yes'));
% 
% cfg              = [];
% cfg.output       = 'pow';
% cfg.channel      = 'MEG';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 1 Hz 
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
% cfg.toi          = 0:0.5:120;                      % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% TFRhann = ft_freqanalysis(cfg, eyesClosedData);
% A.TFRhann = TFRhann; 
% cfg = [];
% % cfg.baseline     = [-0.5 -0.1];
% % cfg.baselinetype = 'absolute';
% % cfg.xlim         = [0.9 1.3];
% cfg.zlim         = [0 6e-26];
% % cfg.ylim         = [A.tfAlpha-2 A.tfAlpha+2];
% cfg.ylim         = [window(1) window(end)];
% cfg.gridscale = 300;
% % cfg.style = 'straight';
% cfg.marker       = 'numbers';
% cfg.channel         = p.megChannels; 
% cfg.colorbar        = 'yes';
% figure
% ft_topoplotTFR(cfg, TFRhann);
% title(sprintf('%s alpha power, 8-14 Hz',und2space(sessionDir)))
% % title(sprintf('%s alpha power, %d ± 2Hz',und2space(sessionDir),A.tfAlpha))

%% topo of tf within alpha range variable color bar 
vals = A.chTfPows(:,window); % isolate tfpows in alpha window 
vals = nanmean(vals,2); % average across alpha window 
figure
meg_topoplot(vals)
colorbar
title(sprintf('%s alpha power, 8-14 Hz',und2space(sessionDir)))

%% topo of tf within alpha range, fixed color bar
vals = A.chTfPows(:,window); % isolate tfpows in alpha window 
vals = nanmean(vals,2); % average across alpha window 
figure
meg_topoplot(vals)
caxis([0 80])
colorbar
title(sprintf('%s alpha power, 8-14 Hz',und2space(sessionDir)))

%% sort channels by tfPow
vals = [];
vals = A.tfPows; % ch x freq x time 
meanVals = nanmean(vals(:,A.tfAlpha,:),3); 
[valsSort,idxSort] = sort(meanVals,'descend'); 
A.channelsRankedAlpha = idxSort; 
A.alphaRanked = valsSort; 

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
figNames = {'tf_EyesClosed','tf_Alpha','tf_Alpha_singleCh','fftAlpha_EyesClosed','fftAlpha_EyesClosed_singleCh','fftApha_Histogram','topoAlphaPower_singleChannel','topoAlphaPower_variableCAxis','topoAlphaPower_fixedCAxis'}; 

end


