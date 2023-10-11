function [E, yhat1, yhat2] = meg_objectiveFunction2(x,data1,data2,t,Fs,paramNames,fitType)
% Define objective function to minimize for fconmin
% August 2023 KT
% Inputs:
%   data: t x 
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

if nargin<5 % if no sampling frequency provided, default to 1000 Hz
    Fs = 1000;
end

if nargin<6
    paramNames = {'intercept1','slope1','amplitude1','phase1',...
              'intercept2','slope2','amplitude2','phase2',...
              'freq'};
end

if nargin<7
    fitType = 'linear';
end

%% Translate starting coefficients into meaningful vars 

% --- Intercept --- 
idx = find(contains(paramNames,'intercept1')); 
intercept1 = x(idx);
idx = find(contains(paramNames,'intercept2')); 
intercept2 = x(idx);

% --- Slope --- 
idx = find(contains(paramNames,'slope1')); 
slope1 = x(idx);
idx = find(contains(paramNames,'slope2')); 
slope2 = x(idx);

switch fitType
    case 'linear2Hz'
        % --- Amplitude ---
        idx = find(contains(paramNames,'amplitude1'));
        amplitude1 = x(idx);
        idx = find(contains(paramNames,'amplitude2'));
        amplitude2 = x(idx);

        % --- Phase ---
        idx = find(contains(paramNames,'phase1'));
        phase1 = x(idx); % in rads??
        idx = find(contains(paramNames,'phase2'));
        phase2 = x(idx);

        % --- Freq ---
        idx = find(contains(paramNames,'freq'));
        freq = x(idx);
        % *pi/(Fs/2);
        % freq = 2; 
end

%% Predicted y 
switch fitType
    case 'linear2Hz'
        phase1_t = rad2t(phase1,freq,Fs);
        phase2_t = rad2t(phase2,freq,Fs); 
        yhat1 = ((slope1/1000 *t) + intercept1)  + ( amplitude1 * sin( (freq*pi/(Fs/2)) * (t + phase1_t * 100 )) );
        yhat2 = ((slope2/1000 *t) + intercept2)  + ( amplitude2 * sin( (freq*pi/(Fs/2)) * (t + phase2_t * 100 )) );
    case 'linear'
        yhat1 = (slope1/1000 * t) + intercept1;
        yhat2 = (slope2/1000 * t) + intercept2;
end

y1 = data1;
y2 = data2;

%% Error to be minimized 
% sum of squared errors 
E = sum ( (yhat1-y1).^2 + (yhat2-y2).^2 ); % CURRENT

%% Normalize (data or errors?)? 
% yhat1_norm = normalize(yhat1); 
% yhat2_norm = normalize(yhat2); 
% y1_norm = normalize(y1); 
% y2_norm = normalize(y2); 

% y1_error = rescale( (yhat1-y1).^2, 0,1); 
% y2_error = rescale( (yhat2-y2).^2, 0,1); 
% 
% figure
% subplot 211
% hold on 
% figureStyle
% histogram((yhat1-y1).^2,100)
% histogram((yhat2-y2).^2,100)
% xlabel('Error')
% ylabel('Count')
% 
% subplot 212
% hold on 
% figureStyle
% bar([1 2 3], [sum((yhat1-y1).^2) sum((yhat2-y2).^2) sum( ((yhat1-y1).^2) + (yhat2-y2).^2)])
% 
% E = sum(y1_error+y2_error);

% E = sum ( (yhat1_norm-y1_norm).^2 + (yhat2_norm-y2_norm).^2 ); 

%% --- OLD ---
% E = sum(rmse(yhat1,y1)) + sum(rmse(yhat2,y2));

%% Calculate rmse 
% diff = sum((y-yhat).^2); 
% E = sqrt(diff/length(y)); 

% Try additive function, no doesn't make sense
% yhat12 = yhat1+yhat2; 
% y12 = y1+y2; 
% E = sum( (yhat12-y12).^2 ); 







