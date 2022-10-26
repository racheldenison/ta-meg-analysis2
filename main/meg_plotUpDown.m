% meg_plotUpDown
% June 29, 2020

% load group channels 
%  load('groupC_TANoise_20Hz.mat')

nChOI = 5; 
saveFigs = 1; 
fieldNames = {'cueT1','cueT2'}; 
p = meg_params('TANoise_Analysis'); 
t = p.tstart:p.tstop;

[sessionNames,subjectNames,ITPCsubject,ITPCsession]  = meg_sessions('TANoise'); 

%% high pass filter

% high pass filter options
Fsample = 1000;
Fhp = 2; % 1, 0.1 high pass frequency
N = []; %16500; % 8250 % filter order, auto calculate if unspecified 
type = 'firws';
direc = 'onepass-zerophase';

for i = 1:numel(sessionNames)
    for f = 1:numel(fieldNames)
        vals = []; 
        vals = groupDavg.(fieldNames{f})(:,:,i); 
        vals = ft_preproc_highpassfilter(vals', Fsample, Fhp, N, type, direc);
        groupDFilter.(fieldNames{f})(:,:,i) = vals';
    end
end

%% get data from top channels 
for i = 1:numel(sessionNames)
    channelsRanked(i,:) = groupC(i).channelsRanked(1:nChOI); 
    channelDir(i,:) = groupC(i).channelDirection(channelsRanked(i,:)); 
end

% top channels, flip direction of channel amp
for i = 1:numel(sessionNames)
    for c = 1:nChOI
        chDir = channelDir(i,c);
        vals_cT1(:,c,i) = groupDavg.cueT1(:,channelsRanked(i,c))*chDir; % time x ch x session
        vals_cT2(:,c,i) = groupDavg.cueT2(:,channelsRanked(i,c))*chDir;
    end
end

%% average and plot 
for i = 1:numel(sessionNames)
    c = 1:nChOI; 
    % for c = 1:nChOI
        figure
        set(gcf, 'Position',  [100, 100, 800, 400])
        hold on
        
        plot(t,vals_cT1(:,c,i),'Color',p.cueColors(1,:),'LineWidth',1.5)
        plot(t,vals_cT2(:,c,i),'Color',p.cueColors(2,:),'LineWidth',1.5)
        %     plot(t,nanmean(vals_cT1(:,1:nChOI,i),2),'Color',p.cueColors(1,:),'LineWidth',2)
        %     plot(t,nanmean(vals_cT2(:,1:nChOI,i),2),'Color',p.cueColors(2,:),'LineWidth',2)
        set(gca,'TickDir','out');
        ax = gca;
        ax.LineWidth = 1.5;
        ax.XColor = 'black';
        ax.YColor = 'black';
        ax.FontSize = 12;
        ylabel('Amplitude (T)')
        xlabel('Time (ms)')
        ax.XLim = [t(1),t(end)]; % trial time
        % ax.XLim = [p.eventTimes(2)-100,p.eventTimes(3)+100]; % targets time
        vline(p.eventTimes,':k',p.eventNames,[-0.01 0.04])
        hline(0,':k')
        title(sprintf('%s channel %d:%d',und2space(sessionNames{i}),c(1),c(end)))
        
        if saveFigs
            filename = sprintf('%s_channel%d_%d.png',sessionNames{i},c(1),c(end));
            foldername = '/Users/kantian/Dropbox/Data/TANoise/MEG/Group/figures/ebi_cue/';
            % print(gcf,'-depsc','-painters',[foldername filename])
        print(gcf,'-dpng',[foldername filename])
        end
        close all
    % end
end

%% ITPC and phase angle 

% groupDavg data: time x ch x trial x session
% fieldtrip input: ch x sample 

% taper     = 'hanning';
foi      = 20;
% t_ftimwin   = 10 ./ foi;
% toi      = tstart/1000:0.01:tstop/1000;
toi = p.eventTimes(2) + abs(p.tstart); 
width = 8;

% field trip trial structure 
for iTrial = 1:192 
    vals = []; 
    vals = groupDavg.cueT1(:,channelsTop,iTrial,i); 
    data.trial{iTrial} = vals'; 
end

for i = 1 % session 
    channelsTop = groupC(i).channelsRanked(1:nChOI); 
    data = squeeze(groupDavg.cueT1(:,channelsTop,:,i)); 
    data = squeeze(nanmean(data,2)); % average channels
    data = permute(data,[2,1,3]); 
    
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
    
    spec = squeeze(spectrum);
    specAmp = abs(squeeze(spectrum));
    itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials

end
