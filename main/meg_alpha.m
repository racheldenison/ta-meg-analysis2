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
foi            = 1:100;
t_ftimwin      = 10 ./ foi;
toi            = t(1)/p.fSample:0.01:t(end)/p.fSample;
tfAmps = [];
window = 5:25; % Hz, to look for peak alpha 

% figures
ytick = 10:10:numel(foi);
xtick = 1:50:numel(toi);

ylims = [min(foi),max(foi)]; 
xlims = [size(toi,1),size(toi,2)]; 

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

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));
figNames = {'tf_EyesClosed','tf_Alpha','fftAlpha_EyesClosed'}; 

end


