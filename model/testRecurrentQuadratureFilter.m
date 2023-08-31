function testRecurrentQuadratureFilter

fs = 1000;                   % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 2.3;              % seconds
t = (0:dt:StopTime-dt)';     % seconds
F = 20;                      % Sine wave frequency (hertz)
amplitude = 0.5; 
verticalShift = 0.5; 

x = amplitude * sin(2*pi*F*t) + verticalShift;  % contrast response (-1 to 1) 
noise = (rand(size(x))-0.5)*0.25; % generate noise 
x = x + noise; % add noise to signal 
targetDur = 100; 
targetTimes = [1050 1350];   % T1 and T2 times 
targetResp = [targetTimes(1):targetTimes(1)+targetDur-1 targetTimes(2):targetTimes(2)+targetDur-1];    
x(targetResp) = 3; 

figure
hold on 
plot(t*1000,x); 
meg_figureStyle

xline(targetTimes(1), 'T1')
xline(targetTimes(2), 'T2')

xlabel('Time (ms)')
ylabel('Contrast (?)') 


%% Quadrature filter 
omega = 20/1000; % preferred temporal frequency, cycles/msec 
lambda = 0.04; % bandwidth 
tau = 1; % bandwidth 

deltaT = 0.1; % msec
RIC = 'r'; 
% rrrr';

n = 1; % number of filters in the cascade
y = recurrentQuadratureFilter_kt(x,lambda,omega,tau,n,RIC,deltaT); 

%% 
for iTrial = 1:5
    % Make signal
    fs = 1000;                   % Sampling frequency (samples per second)
    dt = 1/fs;                   % seconds per sample
    StopTime = 2.3;              % seconds
    t = (0:dt:StopTime-dt)';     % seconds
    F = 20;                      % Sine wave frequency (hertz)
    amplitude = 0.5;
    verticalShift = 0.5;

    x = amplitude * sin(2*pi*F*t) + verticalShift;  % contrast response (-1 to 1)
    noise = (rand(size(x))-0.5)*1; % generate noise
    x = x + noise; % add noise to signal
    targetDur = 100;
    targetTimes = [1050 1350];   % T1 and T2 times
    targetResp = [targetTimes(1):targetTimes(1)+targetDur-1 targetTimes(2):targetTimes(2)+targetDur-1];
    x(targetResp) = 3;

    n = 5; % number of filters in the cascade
    figure(iTrial)
    subplot 211
    RIC = repmat('r',1,n);
    y = recurrentQuadratureFilter_kt(x,lambda,omega,tau,n,RIC,deltaT);
    plot(squeeze(y),'LineWidth',1.5)
    ylabel('Real')

    subplot 212
    RIC = repmat('i',1,n);
    y = recurrentQuadratureFilter_kt(x,lambda,omega,tau,n,RIC,deltaT);
    plot(squeeze(y),'LineWidth',1.5)
    ylabel('Imaginary')
    xlabel('Time (ms')
end





