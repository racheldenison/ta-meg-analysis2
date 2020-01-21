function p = meg_params(expt)

% function p = meg_params(expt)
%
% INPUT
%   expt: string expt type ('TA2') 
% OUTPUT
% p
%   parameters
%
% Karen Tian
% January 2020

switch expt
    case 'TA2'
        % timing 
        p.eventTimes = [0 1050 1300 2100]; % accounting for 50ms photodiode delay
        p.eventNames = {'precue','T1','T2','response cue'};
        p.tstart = -500; % -1000; 
        p.tstop = 2800; % 2300;   
        p.prestim = 0.5; 
        p.poststim = 2.8; 
        
        % channels
        p.megChannels = 1:157; 
        
    otherwise 
        disp('expt type not recognized')
end