function [channelsRanked, C] = meg_loadChannels(matDir, channelSelectionType)

% function [channelsRanked, C] = meg_loadChannels(matDir, channelSelectionType)
%
% INPUTS
% matDir
%   string, full path to mat directory
%
% channelSelectionType
%   method used to sort the channels, corresponds to an existing file name, 
%   'channels_[type name]', e.g. 'channels_peakprom'
%   availabble: 'peakprom', 'classweights', '20Hz_ebi'
%
% OUTPUTS
% channelsRanked 
%   vector, channel numbers rank ordered by the sorting metric
%
% C
%   channel structure returned by channel sorting function. see
%   meg_sortChannels.m
%
%
% April 2020

switch channelSelectionType
    case 'Pk_avgProm'
        %%% update to give channelsRanked and C outputs
        chFile = sprintf('%s/Pk_avgProm.mat',matDir);
        load(chFile)
        C = Pk; 
        channelsRanked = C.idxDirProm; 
        % selectedChannels = Pk.passCh';
        % channelDir = Pk.promDir(selectedChannels); % positive or negative peak
    case {'peakprom','classweights'}
        chFile = sprintf('%s/channels_%s.mat',matDir,channelSelectionType);
        load(chFile)
        channelsRanked = C.channelsRanked;
    case '20Hz_ebi'
        chFile = sprintf('%s/channels_%s.mat',matDir,channelSelectionType);
        C = load(chFile);
        channelsRanked = C.channelsRanked;
    otherwise
        error('channelSelectionType not recognized')
end

