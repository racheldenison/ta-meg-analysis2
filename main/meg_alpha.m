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
p = meg_params('TA2'); % get parameters 
data = data(:,p.megChannels); % read only meg data

%% params 

% timing 
sz = size(data);
t = 1:sz(1); 

taper          = 'hanning';
foi            = 1:25;
t_ftimwin      = 10 ./ foi;
toi            = t(1)/p.fSample:0.01:t(end)/p.fSample;
tfAmps = [];
window = 5:25; % Hz, to look for peak alpha 

% figures
ytick = 10:10:numel(foi);
xtick = 1:50:numel(toi);

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 

load('data_hdr.mat')

nSamples = t(end); 
nfft = 2^nextpow2(nSamples);

%% time freq analysis

y = fft(data,nfft)/nSamples; % freq x trial
f = p.fSample/2*linspace(0,1,nfft/2+1); % frequencies
amps = 2*abs(y(1:nfft/2+1,:,:)); % amp of freq x trial
ampsMean = nanmean(amps,3);

% fieldtrip wavelet convolution
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
meanTfAmps = squeeze(nanmean(tfAmps,1));

% store variables
A.y = y;
A.f = f;
A.amps = amps;
A.ampsMean = ampsMean;
A.tfAmps = tfAmps;
A.meanTfAmps = meanTfAmps;

%% spectrogram plot 

figure
hold on
imagesc(A.meanTfAmps)
xlim(xlims)
ylim(ylims)
xlabel('time (s)')
ylabel('frequency (Hz)')
colorbar
title(sprintf('Eyes Closed %s ', und2space(sessionDir)))

%% tf output alpha peak 

meanmeanTfAmps = nanmean(A.meanTfAmps,2); 
[peak,peakIdx] = max(meanmeanTfAmps(window)); 
A.tfAlpha = window(peakIdx); 

figure
plot(meanmeanTfAmps)
hold on
xlim([window(1) window(end)])
vline(A.tfAlpha,'k')
title(sprintf('tf alpha %d Hz %s',A.tfAlpha,und2space(sessionDir)))

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
xlabel('frequency (Hz)')
ylabel('amp')

% plot peak frequency within window
A.alpha = A.f(peakIdx)+window(1); 
vline(A.alpha,'k')
title(sprintf('%s alpha %.2f Hz',und2space(sessionDir),A.alpha))

%% field trip multiploter tester

cfg                     = [];
cfg.dataset             = eyesClosedFile;
eyesClosedData          = ft_preprocessing(struct('dataset',eyesClosedFile,...
    'channel','MEG','continuous','yes'));

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:30;                         % analysis 2 to 30 Hz in steps of 1 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = 0:0.5:120;                      % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg, eyesClosedData);
A.TFRhann = TFRhann; 

%% plot topo tfrhann
cfg = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'absolute';
% cfg.xlim         = [0.9 1.3];
% cfg.zlim         = [0 6e-26];
cfg.ylim         = [A.tfAlpha-2 A.tfAlpha+2];
cfg.gridscale = 300;
% cfg.style = 'straight';
cfg.marker       = 'numbers';
cfg.channel         = p.megChannels; 
cfg.colorbar        = 'yes';
figure
ft_topoplotTFR(cfg, TFRhann);
title(sprintf('%s alpha power, %d ± 2Hz',und2space(sessionDir),A.tfAlpha))

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
figNames = {'tf_EyesClosed','tf_Alpha','fftAlpha_EyesClosed','topoalphapower'}; 

end


