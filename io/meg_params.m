function p = meg_params(expt)

% function p = meg_params(expt)
%
% INPUT
%   expt: string expt type ('TA2') 
% OUTPUT
% p
%   parameters, including timing, channel info
%
% Karen Tian
% January 2020

p.fSample = 1000; % frame rate, Hz
p.megChannels = 1:157; % 1 indexed

switch expt
    
    case 'TA2'
        % timing 
        p.eventTimes = [0 1050 1300 2100]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; % -1000; 
        p.tstop = 2800; % 2300;   
        p.prestim = 0.5; 
        p.poststim = 2.8; 
   
    case 'TANoise'
        % add timing
        p.eventTimes = [0 1000 1250 2100];
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; 
        p.tstop = 3100; 
        p.prestim = 0.5; 
        p.poststim = 3.1; 
    otherwise 
        disp('expt type not recognized')
end