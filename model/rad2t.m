function val = rad2t(val,freq,Fs)
% function val = rad2t(val,freq,Fs)
% Conversion from rad space to t space 
% val is vector of values to convert 
% freq is frequency 
% Fs is sampling frequency

%% Check inputs and set defaults 
if nargin<2 
    freq = 2; 
    warning('No frequency specified. Default to 2 Hz.')
end
if nargin<3 
    Fs = 1000; 
    warning('No sampling frequency specified. Default to 1000 Hz.')
end

%% Convert vals from rad space to t space 
val = (val*Fs) / (2*freq*pi*100); 


