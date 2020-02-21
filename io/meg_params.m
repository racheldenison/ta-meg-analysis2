function p = meg_params(expt)

% function p = meg_params(expt)
%
% INPUT
%   expt: string expt type ('TA2_Preproc','TA2_Analysis','TANoise_Preproc','TANoise_Analysis') 
% OUTPUT
% p
%   parameters, including timing, channel info
%
% Karen Tian
% January 2020

p.fSample = 1000; % frame rate, Hz
p.megChannels = 1:157; % 
p.eventNames = {'precue','T1','T2','response cue'};

switch expt
    
    case 'TA2_Preproc'
        % timing 
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.tstart = -200; % -1000; 
        p.tstop = 2300; % 2300;   
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        
    case 'TA2_Analysis'
        % timing 
        p.eventTimes = [0 1050 1350 2300]; % accounting for 50ms photodiode delay
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
   
    case 'TANoise_Preproc'
        % add timing
        p.eventTimes = [0 1000 1250 2300];
        p.tstart = -200; 
        p.tstop = 2300; 
        p.prestim = 0.2; 
        p.poststim = 2.3; 
        
    case 'TANoise_Analysis'
        % add timing
        p.eventTimes = [0 1000 1250 2300];
        p.tstart = -500; 
        p.tstop = 3900; 
        p.prestim = 0.5; 
        p.poststim = 3.9; 
    otherwise 
        disp('expt type not recognized')
end