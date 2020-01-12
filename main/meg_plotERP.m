function fH = meg_plotERP(data,selectedChannels,plotSingleTrial,plotAvgTrial,plotAvgChannel)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% data
%   structure of condition cells containing
%       data matrix, time x channels x trials
% selectedChannels
%   vector, default 1:157 meg channels
% plotSingleTrial by channel by condition
%   1 or 0, default 1
% plotAvgTrial
%   1 or 0, default 1
% plotAvgChannel
%   1 or 0, default 1
%
% Karen Tian
% January 2020

%% args

if nargin<5
    plotAvgChannel = 1;
end
if nargin<4 % default plots avg trial
    plotAvgTrial = 1;
end
if nargin<3 % default plots single trial
    plotSingleTrial = 1;
end
if nargin<2 % default selects channels 1:3
    selectedChannels = 1:3;
    disp('selectedChannels not specified, default 1:3')
end
if nargin<1
    load('data.mat'); % load dummy data
end

%% setup
saveFigs = 0;

fieldName = fieldnames(data);
nFields = numel(fieldName);
nChannels = numel(selectedChannels);

eventTimes = [0 1000 1250 2100]+500; % check timing
eventNames = {'precue','T1','T2','response cue'};

t = 1:size(data.(fieldName{1}),1);
xlims = [min(t),max(t)];

%% checks

% check data
sz = size(data.(fieldName{1}));
if numel(sz)~=3
    error('data is expected to be 3-dimensional, with trials as the last dimension')
end
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% plots ERP: fig per channel, subplot per condition, single trials 

if plotSingleTrial
    
    for iC = selectedChannels
        figure
        for iF=1:nFields
            vals = data.(fieldName{iF});
            subplot (nFields,1,iF)
            hold on
            plot(t, squeeze(vals(:,iC,:)))
            title(sprintf('%s',fieldName{iF}))
            xlim(xlims)
            xlabel('time (ms)')
            ylabel('amplitude')
            vline(eventTimes,'k',eventNames)
        end
        rd_supertitle2(sprintf(sprintf('channel %d',iC)))
    end
    
end

%% plot ERP average trial by channel

if plotAvgTrial
    
    figure
    for iC=1:nChannels
        for iF=1:nFields
            vals = data.(fieldName{iF});
            meanTrial = nanmean(vals(:,iC,:),3);
            meanChannel = nanmean(meanTrial,2);
            subplot (nChannels,1,iC)
            hold on
            plot(t, meanChannel)
        end
        xlim(xlims)
        vline(eventTimes,'k',eventNames)
        xlabel('time (ms)')
        ylabel('amplitude')
        title(sprintf('channel %d',iC))
    end
    legend(fieldName)
    rd_supertitle2('avg ERP')
    
end


%% plot ERP average channel

if plotAvgChannel
    
    figure
    for iF=1:nFields
        vals = data.(fieldName{iF});
        meanTrial = nanmean(vals(:,selectedChannels,:),3);
        meanChannel = nanmean(meanTrial,2);
        hold on
        plot(t, meanChannel)  
    end
    xlim(xlims)
    vline(eventTimes,'k',eventNames)
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('channel %d ',selectedChannels))
    legend(fieldName)
    rd_supertitle2('avg ERP')
    
end


%% return figure handle



