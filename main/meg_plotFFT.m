function [A,fH,figNames] = meg_plotFFT(data,p,selectedChannels,selectedFreq)

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

p = meg_params('TANoise_ITPCsession8');

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
nfft = 2^nextpow2(nSamples);
Fsample = 1000; 

cueColors = [122/255 142/255 194/255; 225/255 124/255 96/255; 128/255 128/255 128/255];  % cueT1, cueT2, neutral 
colorAlpha = 0.75; % transparency for plots 

%% generate noisy sine wave 
x = 0:.1:5000; 
y = sin(x); 
idx = 1:1000;
noiseAmplitude = 1; 
noisy_y = y + noiseAmplitude * randn(1, length(y)); % faux noisy sine wave

figure
plot(x(idx), noisy_y(idx))

%% setup fft
% Fs = 1000; 
% freqRes = 0.1; % Fs/nfft; 

nfft = 20; % 2^nextpow2((Fs/freqRes)-1);

%% do fft
data = noisy_y(idx); 
nSamples = length(data); 


y = fft(data,nfft)/nSamples; % freq x trial
f = Fs/2*linspace(0,1,nfft/2+1); % frequencies
amps = 2*abs(y(1:nfft/2+1)); % amp of freq x trial
ampsMean = nanmean(amps,3);

figure
plot(f,amps,'-o','LineWidth',2,'Color',p.cueColors(1,:))
title(sprintf('nfft = %d, nSamples = %d',nfft,length(idx)))
titleText = sprintf('FFT_nfft%d_nSamples%d',nfft,length(idx)); 

% plot check 
set(gca,'TickDir','out');
ax = gca;
box off
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
ylabel('Amplitude')
xlabel('Frequency')

print(sprintf('%s.png',titleText),'-dpng')

%% 

freqRes = Fs/nfft; 

fres = 25;      %desired frequency resolution
N = length(U);  %length of the signal
nfft = 2^nextpow2((fs/fres)-1); %number of spectral lines --> power of 2 because radix-2-algorithm                     
fRes = fs/(nfft+1); %real frequency resolution depending on nfft
ns = floor(N/nfft); %number of sections
%Frequency axis
fax = linspace(0,fs,nfft+1);
fax = fax(1:nfft/2+1);
%deviding the data into ns-sections with length of nfft
U = reshape(U(1,1:ns*nfft), [nfft, ns])';
%window function
win = hann(nfft,'Periodic')';   %creating a hanning window with nfft values
W = repmat(win,ns,1);       %creating "window" matrix to multiply with each section
k=nfft/sum(win,2);    %correction coefficient of window function (not sure if correct)
Uw = U.*W;          %values of the section multiplied by window function
%Calculation of the frequency spectrum
dim=2;       %defining direction for fft()
UF = fft(Uw,nfft,dim);
%Calculating the amplitude
AU = k*abs(UF);   %Multiplication with correction coefficient for window function
AU = (1/(ns*nfft))*sum(AU);   %averaging the sections
AU = 2*AU(1,1:nfft/2+1);      %multiplying the values by 2 --> two spectral lines make up one amplitude
                              %only the first nfft/2+1 values needed because of multiplication by 2
AU(1) = AU(1)/2;

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