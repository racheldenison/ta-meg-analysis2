function p = meg_params(expt)

% function p = meg_params(expt)
%
% INPUT
%   expt: string expt type ('TA2_Preproc','TA2_Analysis','TANoise_Preproc','TANoise_Analysis','TANoise_ITPC','TANoise_ITPC8','Cupcake') 
% OUTPUT
% p
%   parameters, including channel info, event timing, trial segmenting 
%
% Karen Tian
% January 2020

p.fSample = 1000; % frame rate, Hz
p.megChannels = 1:157; % 
p.channelForSaturatingChannels = 168;

% plot colors
p.cueColors = [122/255 142/255 194/255;  % cueT1 blue
    225/255 124/255 96/255; % cueT2 red
    128/255 128/255 128/255; % cueNeutral grey
    157/255 135/255 212/255]; % difference purple
p.colorAlpha = 0.75; % transparency
p.cueErrorBarColors = [163/255 180/255 216/255; 232/255,168/255,154/255];

switch expt    
    case 'TA2_Preproc'
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -200; % -1000; 
        p.tstop = 2300; % 2300;   
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        
    case 'TA2_Analysis'
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
        p.precueChannel = 168; 
        p.blankChannel = 167;  
        p.cueColors = [122/255 142/255 194/255; 225/255 124/255 96/255; 128/255 128/255 128/255; 157/255 135/255 212/255];  % cueT1, cueT2, neutral, difference
        p.colorAlpha = 0.75; % transparency for plots 
   
    case 'TANoise_Preproc'
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -200; 
        p.tstop = 2300; 
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        
    case 'TANoise_Analysis'
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        
    case 'TANoise_ITPC'
        p.eventTimes = [0 1050 1350 2300];
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -2000;
        p.tstop = 6000;
        p.prestim = 2.0;
        p.poststim = 6.0;
        p.precueChannel = 168;
        p.blankChannel = 167;
        p.trialDefTrig = [p.precueChannel,p.blankChannel];
        p.trialTime = 8000; % ms
        
    case 'TANoise_ITPCsession8'
        p.eventTimes = [0 1050 1350 2300];
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -2000;
        p.tstop = 5000;
        p.prestim = 2.0;
        p.poststim = 5.0;
        p.precueChannel = 168;
        p.blankChannel = 167;
        p.trialDefTrig = [p.precueChannel,p.blankChannel];
        p.trialTime = 7000; % ms
        
    case 'Cupcake'
        p.tstart = -1000; 
        p.tstop = 2000; 
        p.prestim = 1; 
        p.poststim = 2;
        triggerChannels = 160:167;
        p.trialDefTrig = triggerChannels(2)+1;  
 
    otherwise 
        disp('expt type not recognized')
end

% all
p.trialDefTrig = [p.precueChannel p.blankChannel];
p.t = p.tstart:p.tstop;
p.trialTime = numel(p.t); % ms


