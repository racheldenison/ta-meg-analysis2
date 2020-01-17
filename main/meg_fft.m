function fH = meg_fft(data)

% MEG_FFT(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
%
% Karen Tian
% January 2020

%% args
if nargin<1
    load('data.mat'); % load dummy data
end

%% checks

fieldName = fieldnames(data);
nFields = numel(fieldName);

% check data
sz = size(data.(fieldName{1}));
if numel(sz)~=3
    error('data is expected to be 3-dimensional, with trials as the last dimension')
end

%% setup 

Fsample = 1000; 
nSamples = sz(1); 
nfft = 2^nextpow2(nSamples); % 2048

%% concatenate data and fft

for iF=1:nFields
    thisF = fieldName{iF}; 
    vals = data.(thisF);
    nTrials = size(vals,3);
    thisVal = [];
    for iT = 1:nTrials
        thisTrial = vals(:,:,iT);
        thisVal = [thisVal; thisTrial]; 
    end
    concatVals.(thisF) = thisVal; 
    
    toFFT = concatVals.(thisF); 
    Y.(thisF) = fft(toFFT,nfft)/nTrials; % scale by number of samples
    f.(thisF) = Fsample/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
    
    amps.(thisF) = 2*abs(Y.(thisF)(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?
    ampsMean.(thisF) = nanmean(amps.(thisF),3);
    
end

%% plot fft

for iF=1:nFields
    thisF = fieldName{iF}; 
    
    figure
    hold on
    set(gcf, 'Position',  [100, 100, 800, 300])
    
    loglog(f.(thisF), ampsMean.(thisF))
    xlim([f.(thisF)(1) f.(thisF)(end)])
    % ylim([10e-15 10e-10])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    title('FFT of concatenated trialed data')
end

