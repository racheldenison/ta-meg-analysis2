% kt_manuscriptFigs13_ITPCFit_fconmin

% August 2023
% Karen 

%% OLD DELETE !

%% settings
saveFigs = 0; 
user = 'kantian'; % kantian karen 

%% Paths
addpath(genpath(pwd))

% Add circular stats toolbox 
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

%% load ITPC data 
load(sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/unused/groupA_ITPCspectrogram_byAtt.mat',user))
p = meg_params('TANoise_ITPCsession8');

%% Setup data 
% y = A * sin(B(x + C)) + D; 
% amplitude = A 
% period = 2*pi/B 
% frequency = 1/period 
% From FFt we know the frequency to fit is 2 Hz , so B = 2/2*pi
% vertical shift = D 
% horizontal shift = C 

% --- Load ITPC data from above ---
sampling = 1:10:7001;
freq = 20; % frequency of interest 

% --- Fit settings ---
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore; % select timing for fit 

% --- All trials --- 
% select ITPC at 20Hz > average to subjects > do fit 
subVals_all = squeeze(A.all.subject(freq,:,:));  % frequency x time x subjects (from avg sessions)
groupVals_all = squeeze(mean(subVals_all,2)); % average across subjects 
% --- Precue T1 trials --- 
subVals_precueT1 = squeeze(A.cueT1.subject(freq,:,:)); 
groupVals_precueT1 = squeeze(mean(subVals_precueT1,2)); 
% --- Precue T2 trials --- 
subVals_precueT2 = squeeze(A.cueT2.subject(freq,:,:)); 
groupVals_precueT2 = squeeze(mean(subVals_precueT2,2)); 

% --- Prepare data as dataset for model ---
ds = dataset(toi', groupVals(toi)); 
vn = {'t','ITPC'}; 
ds = set(ds,'VarNames',vn);

%% 
