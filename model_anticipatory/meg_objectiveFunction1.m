function [E, yhat] = meg_objectiveFunction1(x,data,t,Fs,paramNames,fitType,freq)
% Define objective function to minimize for fconmin
% August 2023 KT
% Inputs:
%   data: values by t, vector 
%   x: initial coefficients for search
%   t: time vector 
%   Fs: sampling frequency, if not defined default to 1000 Hz 
% Ouputs 
%   E: root mean square error 
%   yhat: predicted model output

% Function: 
% y = A * sin(B(x + C)) + D; 
% amplitude = A 
% period = 2*pi/B 
% frequency = 1/period 
% From FFt we know the frequency to fit is 2 Hz , so B = 2/2*pi
% vertical shift = D 
% horizontal shift = C 

%% 
% For 2 Hz frequency, at 1000 Hz sampling rate, then period = 500 samples
% B = 2*pi/period
% B = 2*pi/(500) = 0.0126

% yhat = sin(t); % basic for testing 
% yhat = x(1) * sin(2*pi/(500)*(t + x(2)*Fs/pi)); 
% period = 2*pi/x(4); % fixed freq to 2 Hz

% fprintf('Fitting %d Hz...',freq)
if nargin<7
    freq = 2; % Hz 
    warning('no frequency specified. Default of 2 Hz...')
end

if nargin<4 % if no sampling frequency provided, default to 1000 Hz
    Fs = 1000;
end

if nargin<5
    paramNames = {'intercept','slope','amplitude','phase'};
    error('Specify parameter names')
end

if nargin<6 % default to linear 
    % fitType = 'linear';
    error('Specify fit type (''linear'' or ''linear2Hz'')')
end

%% Translate x into meaningful vars 
% --- Intercept --- 
idx = find(contains(paramNames,'intercept')); 
intercept = x(idx);

% --- Slope --- 
idx = find(contains(paramNames,'slope')); 
slope = x(idx);

switch fitType
    case 'linear2Hz'
        % --- Amplitude ---
        idx = find(contains(paramNames,'amplitude'));
        amplitude = x(idx);

        % --- Phase ---
        idx = find(contains(paramNames,'phase'));
        phase = x(idx);

        % --- Freq (fixed) ---
        % idx = find(contains(paramNames,'freq'));
        % freq = x(idx);
        % *pi/(Fs/2);
end

% if size(x,2) > 2
%     amplitude = x(3); 
% end
% if size(x,2) > 3
%     phase = x(4)*Fs/pi; % in t space 
% end
% if size(x,2) > 4  
%     freqB = x(5)*pi/500; 
% end

%% Predicted y
% slope/1000 & phase*100 are needed to scale search space to
% approximately the same scale for different parameters 
switch fitType
    case 'linear2Hz'
        phase_t = rad2t(phase,freq,Fs);
        yhat = ((slope/1000 * t) + intercept)  + ( amplitude * sin( (freq*pi/(Fs/2)) * (t + phase_t*100 )) );
    case 'linear'
        yhat = (slope/1000 * t) + intercept;
end

y = data;

% if size(x,2) > 4  
%     yhat = ((slope*t) + intercept)  + ( amplitude * sin( freqB * (t + phase )) );
% else
%     yhat = ((slope*t) + intercept);
% end
% y = data;

%% Calculate rmse 
% diff = sum((y-yhat).^2); 
% E = sqrt(diff/length(y)); 

% E = sum(rmse(yhat,y));
% E = rmse(yhat,y,'all');

E = sum ( (yhat-y).^2 ); % divided by samples? 







