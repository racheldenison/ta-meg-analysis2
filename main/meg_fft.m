function [f P1] = meg_fft(data)
% function meg_fft(data)
% 
% data should be t x 1 vector 
% October 2023 

%% fft setup 
% Need at least 2 points to estimate a sin
% If N is the number of points in a signal 
% N/2 is the maximum frequency that can be sampled (i.e. Nyquist frequency)

plotFigs = 0; 

% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

%% Check data
% data = data2.cueT1.session(tIdx,1); 
% 
% x = 1:971; freq = 2; 
% data = sin((freq*pi/(Fs/2)) * x);

%% Zero pad input data 
dataOriginal = data; 
padSize = 0; % 1024*2-numel(dataOriginal); 
dataPadded = padarray(dataOriginal,[padSize 1],0,'post');

%% Setup 
% Fmax = 1/(2*dt); % Nyquist frequency 
Fs = 1000;                      % Sampling frequency                    
T = 1/Fs;                       % Sampling period       
L = length(dataPadded);         % Length of signal
y = fft(dataPadded);            % FFT the signal -- 2-sided frequency, complex
P2 = abs(y/L);                  % 2-sided PSD -- 2-sided frequency, real magnitude
P1 = P2(1:L/2+1);               % 1-sided PSD -- 1-sided frequency, real magnitude
P1(2:end-1) = 2*P1(2:end-1);    % Power, scaled to match time series magnitude (*)
df = Fs/length(dataPadded);     % frequency resolution 
f = Fs*(0:(L/2))/L;             % frequencies 

%% Plot fft 
if plotFigs
    figure
    subplot 121
    hold on
    meg_figureStyle
    plot(f,P1,'-x','LineWidth',2,'Color',p.cueColors(1,:))
    ylabel('Amplitude')
    xlabel('Frequency (Hz)')
    xlim([0 5])

    subplot 122
    hold on
    meg_figureStyle
    plot(log(f),log(P1),'-x','LineWidth',2,'Color',p.cueColors(1,:))
    ylabel('Log(Amplitude)')
    xlabel('Log(Frequency) (Hz)')
    xlim([0 5])
end

%% Plot 1/f 
% lightBLUE = [189 228 255]/255;
% darkBLUE  = [24 38 186]/255;
% blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
% 
% x = 0.1:0.1:10; 
% 
% figure
% subplot 121 
% hold on 
% meg_figureStyle
% baseline = 0; 
% exponent = linspace(0.1,1,7);
% N = 7; % # of curves to be plotted 
% for i = 1:N
%     val_1f{i} = baseline + 1./(x.^exponent(i)); % 1/f data 
% end
% for i = 1:N
%     plot(x,val_1f{i},'color',blueGRADIENTflexible(i,N),'LineWidth',1); %loop over all curves in plot
% end
% ylabel('log(Amplitude)')
% xlabel('log(Frequency) (Hz)')
% xlim([0 5])
% legend( strsplit(num2str(exponent)) )
% title('Exponent')
% 
% subplot 122 
% hold on 
% meg_figureStyle
% baseline = linspace(0.1,1,7);
% exponent = 1;
% N = 7; % # of curves to be plotted 
% for i = 1:N
%     val_1f{i} = baseline(i) + 1./(x.^exponent); % 1/f data 
% end
% for i = 1:N
%     plot(x,val_1f{i},'color',blueGRADIENTflexible(i,N),'LineWidth',1); %loop over all curves in plot
% end
% ylabel('Amplitude')
% xlabel('Frequency (Hz)')
% xlim([0 5])
% legend( strsplit(num2str(baseline)) )
% title('Baseline')


