function val = t2rad(val,freq,Fs)
% function val = t2rad(val,freq,Fs)
% Conversion from t space to rad space 
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

%% Convert vals from t space to rad space 
val = (val*100)*(2*freq*pi)/Fs;
% or 4 times?? 



