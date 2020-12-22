function [amps,freqs] = meg_FFTparams(data,Fs,df)

% MEG_PLOTTF(data,p,selectedChannels,selectedFreq)
% computes fft based on desired frequency resolution and sampling rate 
%
% INPUT
% data
%   structure of condition cells containing
%   data matrix, time x channels x trials
% Fsamplex
%   frames per second (e.g. 1000 Hz)
% df
%   desired frequency resolution (e.g. 1 Hz)
% selectedFreq  
%   vector of frequencies of interest
% OUPUT
% frequencies
%   0:df:df*(N/2)
% fH
%   figure handle
% figNames 
%   figure names 
%
% Karen Tian
% December 2020 

nSamples = size(data,numel(size(data)));
nfft = Fs/df; 

% nfft = 100; % 2^nextpow2(nSamples);
% f = Fs/2*linspace(0,1,nfft/2+1); % frequencies

freqs = 0:df:df*(nfft/2); 
    
y = fft(data,nfft)/nSamples; % freq x trial
amps = 2*abs(y(1:nfft/2+1)); % amp of freq


fH = sort(double(findobj(0,'Type','figure')));
% figNames = {'TF','TFnorm','FFTAllCh','FFTSelectCh','Freq','TF_diffT2T1'}; 
figNames = {'TF','TFnorm','FFTAllCh','FFTSelectCh','Freq'}; 



