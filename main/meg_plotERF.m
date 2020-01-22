function [fH,figNames] = meg_plotERF(data,selectedChannels,plotSingleTrial,plotAvgTrial,plotAvgChannel)

% MEG_PLOTERF(data,selectedChannels,plotSingleTrial,plotAvgTrial)
%
% INPUTS
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
% OUPUT
%   fH  
%       figure handle
%   figNames
%       figure names
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
    selectedChannels = [20,23,36,43,60];
    disp('selectedChannels not specified, default [20,23,36,43,60]')
end
if nargin<1
    load('D3.mat'); % load dummy data
    data = D3; 
end

%% setup

fieldName = fieldnames(data);
nFields = numel(fieldName);
nChannels = numel(selectedChannels);

p = meg_params('TA2');
t = p.tstart:p.tstop;
xlims = [min(t),max(t)];

%% checks

% check data
sz = size(data.(fieldName{1}));
nTrials = sz(3); 
nAllChannels = sz(2); 
if numel(sz)~=3
    error('data is expected to be 3-dimensional, with trials as the last dimension')
end
    
% check channels 
if ~isnumeric(selectedChannels)
    error('selectedChannels is expected to be a vector of channel numbers')
end

%% plots ERF: fig per channel, subplot per condition, single trials 

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
            vline(p.eventTimes,'k',p.eventNames)
        end
        rd_supertitle2(sprintf(sprintf('channel %d',iC)))
    end
    
    % test peaks
%     spacer = 40;
%     figure
%     for iF=1:nFields
%         vals = data.(fieldName{iF});
%         subplot (nFields,1,iF)
%         hold on
%         for iC = 1:nAllChannels
%             meanTrial = nanmean(vals(:,iC,:),3);
%             plot(t, abs(meanTrial) + spacer*iC)
%         end
%         title(sprintf('%s',fieldName{iF}))
%         xlim(xlims)
%         xlabel('time (ms)')
%         ylabel('amplitude')
%         vline(p.eventTimes,'k',p.eventNames)
%     end
    
end

%% single channel find peaks 

% find peaks
% figure 
% for iF=1:nFields
%     vals = data.(fieldName{iF});
%     subplot (nFields,1,iF)
%     hold on
%     spacer = 200; 
%     for iC = 1:20
%         meanTrial = abs(nanmean(vals(:,iC,:),3));
%         [pks,locs,w,prom] = findpeaks(meanTrial,t); 
%         pkinfo = [pks locs' w' prom];
%         [pSort,pIdx] = sort(prom,'descend');
%         
%         plot(t, meanTrial + iC*spacer)
%         topPeaks = pIdx(1:3); 
%         plot(locs(topPeaks),pks(topPeaks)+iC*spacer,'.','MarkerSize',30,'Color',[0,1,0])
%         
%     end
%     title(sprintf('%s',fieldName{iF}))
%     xlim(xlims)
%     xlabel('time (ms)')
%     ylabel('amplitude')
%     vline(p.eventTimes,'k',p.eventNames)
% end

%% plot ERF average trial by channel

if plotAvgTrial
    
    figure
    for iC=selectedChannels
        for iF=1:nFields
            vals = data.(fieldName{iF});
            meanTrial = nanmean(vals(:,iC,:),3);
            meanChannel = nanmean(meanTrial,2);
            idx = find(selectedChannels==iC);
            subplot (nChannels,1,idx)
            hold on
            plot(t, meanChannel)
        end
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        xlabel('time (ms)')
        ylabel('amplitude')
        title(sprintf('channel %d',iC))
    end
    legend(fieldName)
    rd_supertitle2('avg ERF')
    
end

%% plot ERF average trial across selected channels

if plotAvgChannel
    
    figure
    set(gcf, 'Position',  [100, 100, 800, 300])
    for iF=1:nFields
        vals = data.(fieldName{iF});
        meanTrial = nanmean(vals(:,selectedChannels,:),3);
        meanChannel = nanmean(meanTrial,2);
        hold on
        plot(t, meanChannel)  
    end
    xlim(xlims)
    vline(p.eventTimes,'k',p.eventNames)
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('channel %d ',selectedChannels))
    legend(fieldName)
    rd_supertitle2('avg ERF')
    
end

%% return figure handle

fH = sort(double(findobj(0,'Type','figure')));

chNames = {}; 
for iC = 1:nChannels
    chNames(iC)= {sprintf('ERF_SingleTrial_Channel%d',selectedChannels(iC))}; 
end
figNames = [chNames,'ERF_avgTrials','ERF_avgTrialsCh'];

