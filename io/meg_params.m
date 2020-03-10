function p = meg_params(expt)

% function p = meg_params(expt)
%
% INPUT
%   expt: string expt type ('TA2_Preproc','TA2_Analysis','TANoise_Preproc','TANoise_Analysis','Cupcake') 
% OUTPUT
% p
%   parameters, including channel info, event timing, trial segmenting 
%
% Karen Tian
% January 2020

p.fSample = 1000; % frame rate, Hz
p.megChannels = 1:157; % 
p.channelForSaturatingChannels = 168;

switch expt
    
    case 'TA2_Preproc'
        % timing 
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -200; % -1000; 
        p.tstop = 2300; % 2300;   
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        p.trialDefTrig = [p.precueChannel,p.blankChannel]; 
        p.trialTime = 2501; % ms 
        
    case 'TA2_Analysis'
        % timing 
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        p.trialDefTrig = [p.precueChannel,p.blankChannel]; 
        p.trialTime = 4401; % ms 
   
    case 'TANoise_Preproc'
        % add timing
        p.eventTimes = [0 1000 1250 2300];
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -200; 
        p.tstop = 2300; 
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        p.trialDefTrig = [p.precueChannel,p.blankChannel]; 
        p.trialTime = 2501; % ms 
        
    case 'TANoise_Analysis'
        % add timing
        p.eventTimes = [0 1000 1250 2300];
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
        p.precueChannel = 168; 
        p.blankChannel = 167; 
        p.trialDefTrig = [p.precueChannel,p.blankChannel]; 
        p.trialTime = 2501; % ms 
        
    case 'Cupcake'
        p.tstart = -1000; 
        p.tstop = 2000; 
        p.prestim = 1; 
        p.poststim = 2;
        triggerChannels = 160:167;
        p.trialDefTrig = triggerChannels(2)+1; 
        p.trialTime = 3001; % ms 
 
    otherwise 
        disp('expt type not recognized')
end

