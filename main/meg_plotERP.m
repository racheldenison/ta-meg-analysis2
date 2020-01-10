function meg_plotERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
% 
% data 
%   structure of condition cells containing
%       data matrix, time x channels x trials
% selectedChannels 
%   vector, default 1:157 meg channels
% plotSingleTrial
%   1 or 0, default 1
% plotAvgTrial
%   1 or 0, default 1
% figDir 
%   path to figure directory 
%
% Karen Tian
% January 2020 

%% args 

% if nargin<5
%     figDir = '/Users/kantian/Dropbox/CarrascoLab/TA2Figs';
% end
if nargin<4 % default plots avg trial 
    plotAvgTrial = 1; 
end
if nargin<3 % default plots single trial 
    plotSingleTrial = 1;  
end
if nargin<2 % default selects channels 1:3 
    selectedChannels = 1:3;    
end
if nargin<1
    load('data.mat'); % load dummy data 
end

saveFigs = 0;
 
fieldName = fieldnames(data);

eventTimes = [0 1000 1250 2100]; 
eventNames = {'precue','T1','T2','response cue'};

%% plot ERP single trial by channel

if plotSingleTrial
    
    for iF=1:numel(fieldName)
        figure
        dataField = data.(fieldName{iF}); 
        t = 1:size(dataField,1); 
        xlims = [min(t),max(t)]; 
        nPlot = 1; 
        if(isnumeric(dataField))
            for iC = selectedChannels
                subplot (length(selectedChannels),1,nPlot)
                hold on
                plot(t, squeeze(dataField(:,iC,:)))
                title(sprintf('%s channel %d',fieldName{iF},iC))
                xlim(xlims)
                vline(eventTimes,'k',eventNames)
                xlabel('time (ms)')
                ylabel('amplitude')
                nPlot = nPlot + 1; 
            end
        end
    end
    
end

%% plot ERP average trial by channel 

if plotAvgTrial
    
    for iF=1:numel(fieldName)
        figure
        dataField = data.(fieldName{iF}); 
        t = 1:size(dataField,1); 
        xlims = [min(t),max(t)]; 
        nPlot = 1; 
        if(isnumeric(dataField))
            for iC = selectedChannels
                subplot (length(selectedChannels),1,nPlot)
                hold on
                meanTrial = nanmean(dataField(:,iC,:),3);               
                plot(t, squeeze(meanTrial))
                title(sprintf('%s channel %d mean',fieldName{iF},iC))
                xlim(xlims)
                vline(eventTimes,'k',eventNames)
                xlabel('time (ms)')
                ylabel('amplitude')
                nPlot = nPlot + 1; 
            end
        end
    end
    
end

%% save figs

if saveFigs
    if ~exist(figDir,'dir')
        mkdir(figDir)
    end
    f = sort(findobj('Type','figure'));
    for iF = 1:numel(f)
        if isnumeric(f)
            figNames{iF} = num2str(f(iF));
        else
            figNames{iF} = num2str(f(iF).Number);
        end
    end
    rd_saveAllFigs(f,figNames,[],figDir);
end


