function y = recurrentQuadratureFilter(x,lambda,omega,tau,n,RIC,deltaT)
% Feb 2023 
% y = recurrentQuadratureFilter(x,lambda,omega,tau,n,RIC,deltaT)
%
% x: 3D input signal with time in the 1st dimension (real or complex
%    valued) 
% y: output signal (real or complex valued: see RIC)
% omega: preferred temporal frequency
% lambda & tau: determine bandwidth
% n: number of filters in the cascade
% RIC: string of 'r', 'i', and 'c' of length n (default: all 'c').
%    Return only the real part or the imaginary part of the fitlered signal
%    from each step of the cascade.
% deltaT: time step (default: 1 msec)
%
% Computes:
%
%    tau dy/dt = -y + lambda (x + i imag(y)) + (1 - lambda) yhat
%    yhat = w y
%    w = 1 + i 2 pi tau omega
%
% Or:
%
%    tau dy_n/dt = -y_n + lambda y_(n-1) + (1 - lambda) yhat
%    yhat = w y_n
%
% Based on Eqs. 51 and 52 in:
% Heeger, D. J. and W. E. Mackey (2018). "ORGaNICs: A Theory of Working
% Memory in Brains and Machines." arXiv preprint arXiv:1803.06288.
% 
% DJH 1/2019

if ~exist('RIC','var')
  RIC = char('c'*ones(1,n));
end
if ~exist('deltaT','var')
  deltaT = 1;
end
debug = 0; 

w = 1 + 1i*2*pi*tau*omega;

T = size(x,1);
y = zeros([T size(x,2) size(x,3) n]);
xn = zeros([T size(x,2) size(x,3)]);
yTmp = zeros([1 size(x,2) size(x,3)]);
yHat = zeros([1 size(x,2) size(x,3)]);
deltaY = zeros([1 size(x,2) size(x,3)]);

xn(:,:,:) = x(:,:,:);
for nn = 1:n
  realX = isreal(xn);
  yn = zeros(size(x));
  for tt = 1:T-1
    yTmp(1,:,:) = yn(tt,:,:);
    yHat(1,:,:) = w * yTmp;
    if realX
      deltaY(:,:,:) = -yTmp + lambda * (xn(tt,:,:) + 1i*imag(yTmp)) + (1-lambda) * yHat;
    else
      deltaY(:,:,:) = -yTmp + lambda * xn(tt,:,:) + (1-lambda) * yHat;
    end
    deltaY(:,:,:) = (deltaT./tau) * deltaY;
    yn(tt+1,:,:) = yTmp + deltaY;
  end
  if strcmp(RIC(nn),'r')
    y(:,:,:,nn) = real(yn);
  elseif strcmp(RIC(nn),'i')
    y(:,:,:,nn) = imag(yn);
  elseif strcmp(RIC(nn),'c')
    y(:,:,:,nn) = yn;
  end
  xn(:,:,:) = y(:,:,:,nn);
end

return

%% Test/debug

if debug 

% clear all; 

% Parameters
deltaT = 0.1; % msec
tau = 1; lambda = 0.04; n = 5; 
omega = 8/1000; % cycles/msec
% RIC = 'ccccc';
RIC = 'rrrrr';
% RIC = 'iiiii';
% RIC = 'icccc';

% Sampling
t = [-100:deltaT:1000-deltaT]';
[xIm,tIm,yIm] = meshgrid([1:3],t,[1:3]);

% Input: sinusoidal flicker
tf = 8/1000; % cycles/msec
flicker = sin(2*pi*tf*tIm);
flicker(:,:,t<0) = 0;
yMaxFlicker = 1.0;

% Input: impulse
impulse = (tIm == 0) / deltaT;
yMaxImpulse = lambda/tau;

% Input: step
step = (tIm >0);
yMaxStep = sqrt(2)/2;

% Pick one
stimulus = impulse; yMax = yMaxImpulse;
% stimulus = flicker; yMax = yMaxFlicker;
% stimulus = step; yMax = yMaxStep;

% Filtered output
y = recurrentQuadratureFilter(stimulus,lambda,omega,tau,n,RIC,deltaT);
yCrop = squeeze(y(:,2,2,:));
yAmp = abs(yCrop);
yIntegral = sum(yCrop) * deltaT;
yAmpIntegral = 0.5 * sum(yAmp) * deltaT;

%% Plot it
figure(1); % clf;
subplot(n+1,1,1);
plot(t,squeeze(stimulus(:,2,2)));
%ylim([0 1/deltaT]);
% drawPublishAxis;
for nn = 1:n
  subplot(n+1,1,nn+1);
  if isreal(yCrop(:,nn))
    plot(t,yCrop(:,nn));
  else
    plot(t,[real(yCrop(:,nn)) imag(yCrop(:,nn))]);
    hold on
    plot(t, yAmp(:,nn), 'k');
    hold off
  end
  %ylim([-yMax yMax]);
  title(['n = ',num2str(nn)]);
  % drawPublishAxis;
end
% print(['Figures/tfilt1.eps'],'-depsc','-painters');

%% Frequency response (this only makes sense for impulse stimulus)
yCrop2 = yCrop(t>=0,:);
F = size(yCrop2,1);
Fs = 1/(deltaT/1000);
frequencies = Fs/2 * linspace(0,1,F/2+1);
yFourierAmp = abs(fft(yCrop2))/F;
yFourierAmp = yFourierAmp(1:F/2+1,:);

% Plot it
figure(2); % clf;
prefFreqs = zeros(1,n);
for nn = 1:n
  subplot(n,1,nn);
  yFourierAmpSc = yFourierAmp(:,nn) / max(yFourierAmp(:,nn));
  prefFreqs(nn) = frequencies(find(yFourierAmpSc == 1));
  loglog(frequencies(2:101), yFourierAmpSc(2:101));
  xlim([1,100])
  ylim([0.01,1]);
  %axis square;
  title(['n = ',num2str(nn)]);
  % drawPublishAxis;
end
prefFreqs
% print(['Figures/tfilt2.eps'],'-depsc','-painters');
end