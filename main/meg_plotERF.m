function [A, fH,figNames] = meg_plotERF(data,p,selectedChannels,plotSingleTrial,plotAvgTrial,plotAvgChannel)

% MEG_PLOTERF(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% INPUTS
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
% p
%   params with expt short type timing info
% selectedChannels
%   vector, default 1:157 meg channels
% plotSingleTrial by channel by condition
%   1 or 0, default 1
% plotAvgTrial
%   1 or 0, default 1
% plotAvgChannel
%   1 or 0, default 1
% 
% OUPUT
%   fH  
%       figure handle
%   figNames
%       figure names
%
% Karen Tian
% January 2020

%% args

if nargin<6
    plotAvgChannel = 1;
end
if nargin<5 % default plots avg trial
    plotAvgTrial = 1;
end
if nargin<4 % default does not plot single trial
    plotSingleTrial = 1;
end
if nargin<3 % default selects channels 1:3
    selectedChannels = [20,23,36,43,60];
    disp('selectedChannels not specified, default [20,23,36,43,60]')
end
if nargin<2
    p = meg_params('TA2_Analysis'); 
    disp('timing parameters not specified, using TA2 analysis')
end
if nargin<1
    load('D3.mat'); % load dummy data
    data = D3; 
end

%% setup

condNames = fieldnames(data);
nConds = numel(condNames);
nChannels = numel(selectedChannels);

t = p.tstart:p.tstop;
xlims = [min(t),max(t)];

cueColors = p.cueColors; 
colorAlpha = p.colorAlpha; 

%% checks

% check data
sz = size(data.(condNames{1}));
if numel(sz)<2 || numel(sz)>3
    error('data is expected to be at least 2-dimensional, with time as the first dimension, and channels second')
elseif numel(sz)==2 
    nTrials = 1; 
elseif numel(sz)==3
    nTrials = sz(3); 
end
nAllChannels = sz(2); 
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% plots ERF: fig per channel, subplot per condition, single trials 

if plotSingleTrial
    
    for iC = selectedChannels
        figure
        for iF=1:nConds
            vals = data.(condNames{iF});
            subplot (nConds,1,iF)
            hold on
            plot(t, squeeze(vals(:,iC,:)))
            title(sprintf('%s',condNames{iF}))
            xlim(xlims)
            xlabel('time (ms)')
            ylabel('amplitude')
            vline(p.eventTimes,'k',p.eventNames)
        end
        rd_supertitle2(sprintf(sprintf('channel %d',iC)))
    end
    
end

%% plot ERF average trial by channel

if plotAvgTrial
    
    figure 
    set(gcf,'Position',[100 100 1500 1200])
    for iC=selectedChannels
        for iF=1:nConds
            vals = data.(condNames{iF});
            meanTrial = nanmean(vals(:,iC,:),3);
            meanChannel = nanmean(meanTrial,2);
            idx = find(selectedChannels==iC);
            subplot (nChannels,1,idx)
            hold on
            figERF = plot(t, meanChannel,'Color',cueColors(iF,:),'LineWidth',2); 
            figERF.Color(4) = colorAlpha; 
        end
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        xlabel('time (ms)')
        ylabel('amplitude')
        title(sprintf('channel %d',iC))
    end
    legend(condNames)
    rd_supertitle2('avg ERF')
    
end

%% plot ERF average trial across selected channels

if plotAvgChannel
    
    figure
    set(gcf, 'Position',  [100, 100, 800, 300])
    for iF=1:nConds
        vals = data.(condNames{iF});
        meanTrial = nanmean(vals(:,selectedChannels,:),3);
        meanChannel = nanmean(meanTrial,2);
        hold on
        figERF = plot(t, meanChannel,'Color',cueColors(iF,:),'LineWidth',1.5);
        figERF.Color(4) = colorAlpha; 
    end
    xlim(xlims)
    vline(p.eventTimes,'k',p.eventNames)
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('channel %d ',selectedChannels))
    legend(condNames)
    rd_supertitle2('avg ERF')
    
end

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));

chNames = {}; 
for iC = 1:nChannels
    chNames(iC)= {sprintf('ERF_SingleTrial_Channel%d',selectedChannels(iC))}; 
end
figNames = [chNames,'ERF_avgTrials','ERF_avgTrialsCh'];

