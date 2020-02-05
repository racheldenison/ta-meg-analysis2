function [fH,figNames] = meg_selectChannels(sessionDir,data)

% MEG_SELECTCHANNELS(sessionDir)
%
% INPUTS
% sessionDir 
% data 
%   matrix time x channel x trial 
% 
% OUPUT
%   C
%       structure of top channels by analysis type 
% or Pk 
% peak analysis 
%
% Karen Tian
% January 2020


%% channel selection method 

selectPeak = 1; % sorts channels by ERF peak
selectThreshold = 1; 
selectAlpha = 0; % sorts channels by alpha power 

%% setup

sessionDir = 'R0817_20181120';
exptShortName = 'TA2';
analStr = 'bietfp';

exptDir = '/Users/kantian/Dropbox/Data/TA2/MEG'; 

p = meg_params('TA2');

fileBase = sessionDirToFileBase(sessionDir, exptShortName);
% 
dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
figDir = sprintf('%s/figures/%s', dataDir);

% behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
% behavFile = dir(sprintf('%s/*.mat', behavDir));
% behav = load(sprintf('%s/%s', behavDir, behavFile.name));

% eyesClosedBase = sessionDirToFileBase(sessionDir, 'EyesClosed');
% eyesClosedFile = sprintf('%s/%s.sqd', dataDir, eyesClosedBase); 
% eyesClosedFileBI = sprintf('%s/%s_bi.sqd', dataDir, eyesClosedBase); 

%% data info

sz = size(data); 
nT = sz(1); 
nCh = sz(2); 
nTrial = sz(3); 

t = p.tstart:p.tstop;
xlims = [min(t),max(t)]; % epoch time in ms

% peak 
vals = nanmean(data,3); % avg across trial 
nPk = 7;

% threshold by max(abs(amp))
toi = p.eventTimes(2):p.eventTimes(3)+300; % times during which to check if max amp crossed  
nChOI = 5; 
percentileCh = (nChOI/nCh)*100; 

%% peak selector

if selectPeak
    % tWindow = p.eventTimes(2):p.eventTimes(3)+300; % window to look for peaks, from T1 to T2+300
    % valsWindow = vals(tWindow,:);
    for iC = 1:nCh
        if nanmean(vals(:,iC),1)==0 % why getting flat channels? check preproc bad channel output
            Pk.peaksPos(iC,:) = NaN;
            Pk.peaksPosVals(iC,:) = NaN;
            Pk.peaksNeg(iC,:) = NaN;
            Pk.peaksNegVals(iC,:) = NaN;
        else
            % positive peaks
            [pks,locs,~,prom] = findpeaks(vals(:,iC));
            [~, idx] = sort(prom,1,'descend'); % idx sorted by desc prominence
            % idx = find(prom>std(prom)*thresh); peak past threshold s
            idx = sort(idx(1:nPk)); % idx ascending of top nPk prom
            peaksPos = t(locs(idx));
            peaksPosVals = pks(idx);
            
            % negative peaks
            [pks,locs,~,prom] = findpeaks(-vals(:,iC));
            [~, idx] = sort(prom,1,'descend');
            idx = sort(idx(1:nPk));
            peaksNeg = t(locs(idx));
            peaksNegVals = -pks(idx);
            
            Pk.peaksPos(iC,:)       = peaksPos;
            Pk.peaksPosVals(iC,:)   = peaksPosVals;
            Pk.peaksNeg(iC,:)       = peaksNeg;
            Pk.peaksNegVals(iC,:)   = peaksNegVals;
        end
    end
    
    % subplot sizing
    nRow = 10;
    nBreak = 35;
    nCol = ceil(nBreak/nRow);
    
    figure
    set(gcf,'Position',[100 100 2400 1200])
    hold on
    nFig = 1;
    for iC = 1:nCh
        if mod(iC,nBreak) == 0 % new figure at breakpoints
            figure
            set(gcf,'Position',[100 100 2400 1200])
            hold on
            nFig = nFig + 1;
        end
        subplot(nRow,nCol,iC-(nBreak*(nFig-1))+1)
        hold on
        plot(t,vals(:,iC),'LineWidth',1)
        plot(Pk.peaksPos(iC,:), Pk.peaksPosVals(iC,:), '.g', 'MarkerSize', 10)
        plot(Pk.peaksNeg(iC,:), Pk.peaksNegVals(iC,:), '.b', 'MarkerSize', 10)
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        yline(0,'--');
        % vline(Pk.peaksPos(iC,:),'g')
        % vline(Pk.peaksNeg(iC,:),'b')
        title(sprintf('channel %d',iC))
        % ylim([min(min(vals)),max(max(vals))])
    end
    
    % store results of peaks analysis
    Pk.t         = t;
    Pk.vals      = vals;
    Pk.nPk       = nPk;
end

%% threshhold selection 

if selectThreshold 
    
    toiY = zeros(length(toi),1); 
    absVals = abs(vals); 
    toiVals = absVals(toi-p.tstart,:); % remember to rereference to 0 for indexing with epoch start time
    toiMax = max(toiVals,[],1); 

    thresholdPrctile = 100-percentileCh; % 96.8153% for 5 channels 
    thresholdVal = prctile(toiMax,thresholdPrctile); % cutoff val
    passCh = find(toiMax > thresholdVal); 
    
    % plot check 
    figure
    plot(toiMax)
    hline(thresholdVal)
    
    % subplot sizing
    nRow = 10;
    nBreak = 35;
    nCol = ceil(nBreak/nRow);
    
    figure
    set(gcf,'Position',[100 100 2400 1200])
    hold on
    nFig = 1;
    for iC = 1:nCh 
        if mod(iC,nBreak) == 0 % new figure at breakpoints
            figure
            set(gcf,'Position',[100 100 2400 1200])
            hold on
            nFig = nFig + 1;
        end
        subplot(nRow,nCol,iC-(nBreak*(nFig-1))+1) % one subplot per channel 
        hold on
        plot(t,absVals(:,iC)) 
        plot(toi,toiY,'k','LineWidth',3) 
        xlim(xlims)
        ylim([0,1.1*max(toiMax)])
        vline(p.eventTimes,'k',p.eventNames)
        yline(0,'--');
        yline(thresholdVal,'r');
        if ismember(iC,passCh)
            yline(thresholdVal,'g'); % if channel passes threshold, turn line green
        end
        title(sprintf('channel %d',iC))
    end
    
    disp(passCh)
end

%%  
%     % try max method 
%     window1 = p.eventTimes(2):p.eventTimes(3); % window to check for first max (T1) 
%     [peak1, idx1] = max(vals(window1,:),[],1); % returns max within defined window, and index relative to window
%     idx1 = idx1 + min(window1) - 1; % convert window to trial time 
% 
%     % T2 mean peak per channel
%     window2 = p.eventTimes(3):p.eventTimes(3)+250; % window to check for second max 
%     [peak2, idx2] = max(vals(window2,:),[],1); % returns max within defined window, and index relative to window
%     idx2 = idx2 + min(window2) - 1; % convert window to trial time 
%     
%     % save 
%     Pk.max1 = peak1; 
%     Pk.idx1 = idx1; 
%     Pk.max2 = peak2; 
%     Pk.idx2 = idx2; 
%     
%     % plot
%     figure
%     set(gcf,'Position',[100 100 2400 1200])
%     hold on
%     nFig = 1;
%     for iC = 1:nCh
%         if mod(iC,nBreak) == 0 % new figure at breakpoints
%             figure
%             set(gcf,'Position',[100 100 2400 1200])
%             hold on
%             nFig = nFig + 1;
%         end
%         subplot(nRow,nCol,iC-(nBreak*(nFig-1))+1)
%         hold on
%         plot(t,vals(:,iC))
%         plot(Pk.peaksPos(iC,:), Pk.peaksPosVals(iC,:), '.g', 'MarkerSize', 10)
%         plot(Pk.peaksNeg(iC,:), Pk.peaksNegVals(iC,:), '.b', 'MarkerSize', 10)
%         xlim(xlims)
%         vline(p.eventTimes,'k',p.eventNames)
%         vline(Pk.idx1(iC),'r')
%         vline(Pk.idx2(iC),'r')
%         hline(0)
%         title(sprintf('channel %d',iC))
%     end
    
%% try a topo 

figure
set(gcf,'Position',[100 100 800 800])
meg_multiplot(vals', [], [], passCh)

%% select based off top alpha channels 

if selectAlpha
    load(sprintf('%s/Alpha.mat',matDir));
    
    power = Alpha.TFRhann.powspctrm;
    powerAvg = nanmean(power,3); % average power across time
    alpha = nanmean(powerAvg(:,Alpha.tfAlpha-2:Alpha.tfAlpha+2),2); % individual alpha ± 2Hz, avg across time
    [C.powSortAlpha,C.chSortAlpha] = sort(alpha,'descend');
    
    figure
    set(gcf,'Position',[100 100 300 500])
    imagesc(alpha); % check
    set(gca,'xtick',[],'ytick',0:10:160)
    ylabel('channel')
    title(sprintf('%s eyes closed alpha %d ± 2Hz \n top 5 channels [%d %d %d %d %d]',und2space(sessionDir),Alpha.tfAlpha,C.chSortAlpha(1:5)))
end

%% save figs

fH = sort(double(findobj(0,'Type','figure')));
% figNames = {'TopChannels_Alpha'};
figNames = {'ERFPeaks1','ERFPeaks2','ERFPeaks3','ERFPeaks4','ERFPeaks5','TopoERF'}
% rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])

%% save top channels 

save(sprintf('%s/Pk.mat',matDir),'Pk')
% save(sprintf('%s/C.mat',matDir),'C')

end



