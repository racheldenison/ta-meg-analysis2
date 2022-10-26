function [fH,figNames,Pk] = meg_channelsPeakSort(sessionDir,data,nChOI,p,passCh)

% MEG_SELECTCHANNELS(sessionDir)
%
% INPUTS
% sessionDir 
% data 
%   matrix time x channel x trial 
% channel selection method 
%   1. ampMax: max(abs(amp)) 
%   2. promAvg: max prominence across two peaks 
%   3. promMax: max prominence averaged from peak in toi1 and peak in toi2
% nChOI 
%   number of channels to select
% 
% Karen Tian
% Feb 2020

%% arguments (add checks) 

selectionMethod = 'promAvg'; % ampMax, promAvg, promMax
Pk.selectionMethod = selectionMethod; 
if nargin<3
    nChOI = 5; 
end

%% setup

% sessionDir = 'R0817_20190625'; % TA2 
% sessionDir = 'R0817_20171212'; % TANoise

exptShortName = 'TANoise';
analStr = 'ebi'; % 'bietfp'

exptDir = sprintf('/Users/kantian/Dropbox/Data/%s/MEG',exptShortName);

fileBase = sessionDirToFileBase(sessionDir, exptShortName);

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);
preprocDir = sprintf('%s/preproc', dataDir);

filename = sprintf('%s/%s_%s.sqd', preprocDir, fileBase, analStr); % *run file* 
dataFile = sprintf('%s/%s_%s_data.mat', matDir, fileBase, analStr);
figDir = sprintf('%s/figures/Channels20Hz', dataDir);

if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% data info

% load(dataFile, 'data')
sz = size(data); 

nT = sz(1); 
nCh = sz(2); 
% nTrial = sz(3); 
% vals = nanmean(data,3); % avg across trial 
vals = data; 
absVals = abs(vals);

% how many channels to select 
percentileCh = (nChOI/nCh)*100; 
thresholdPrctile = 100-percentileCh; % 96.8153% for 5 channels

% p = meg_params(sprintf('%s_Analysis',exptShortName)); 
t = p.tstart:p.tstop;
xlims = [min(t),max(t)]; % epoch time in ms

toi1 = p.eventTimes(2):p.eventTimes(3);
toi2 = p.eventTimes(3):p.eventTimes(3)+length(toi1)-1; 

toiIdx1 = toi1+abs(p.tstart); % reference toi idx to epoch timing
toiIdx2 = toi2+abs(p.tstart); 

toiY1 = zeros(length(toi1),1); % for visualizing toi on x axis 
toiY2 = zeros(length(toi2),1);

% toiVals1 = absVals(toiIdx1,:);
% toiVals2 = absVals(toiIdx2,:);

%% channel selection method: 
%   1. maxAmp: max(abs(amp)
%   2. avgProm: max prominence averaged from two peaks
%   3. maxProm: max prominence across two peaks 
switch selectionMethod
    %% 
    case 'promAvg'
    nPk = 1; % find 1 pos 1 neg peak for T1 window, T2 window 
    for iC = 1:nCh
        if nanmean(vals(:,iC),1)==0 % why getting flat channels? check preproc bad channel output
            Pk.peaksPos1(iC,:) = NaN;
            Pk.peaksPosVals1(iC,:) = NaN;
            Pk.peaksNeg1(iC,:) = NaN;
            Pk.peaksNegVals1(iC,:) = NaN;
            Pk.peaksPos2(iC,:) = NaN;
            Pk.peaksPosVals2(iC,:) = NaN;
            Pk.peaksNeg2(iC,:) = NaN;
            Pk.peaksNegVals2(iC,:) = NaN;
        else
            % positive peak T1
            [pksPos1,locsPos1,~,promsPos1] = findpeaks(vals(toiIdx1,iC));
            [~, idxPos1] = sort(promsPos1,1,'descend'); % idx sorted by desc prominence
            % idx = find(prom>std(prom)*thresh); peak past threshold s
            idxPos1 = sort(idxPos1(1:nPk)); % idx of top nPk prom
            promPos1 = promsPos1(idxPos1); 
            peaksPos1 = toi1(locsPos1(idxPos1)); % idx of peak1
            peaksPosVals1 = pksPos1(idxPos1); % value of peak1
            
            % negative peak T1
            [pksNeg1,locsNeg1,~,promsNeg1] = findpeaks(-vals(toiIdx1,iC));
            [promNeg1, idxNeg1] = sort(promsNeg1,1,'descend');
            idxNeg1 = sort(idxNeg1(1:nPk));
            promNeg1 = promsNeg1(idxNeg1); 
            peaksNeg1 = toi1(locsNeg1(idxNeg1));
            peaksNegVals1 = -pksNeg1(idxNeg1);
            
            % positive peak T2
            [pksPos2,locsPos2,~,promsPos2] = findpeaks(vals(toiIdx2,iC));
            [promPos2, idxPos2] = sort(promsPos2,1,'descend'); % idx sorted by desc prominence
            % idx = find(prom>std(prom)*thresh); peak past threshold s
            idxPos2 = sort(idxPos2(1:nPk)); % idx ascending of top nPk prom
            promPos2 = promsPos2(idxPos2); 
            peaksPos2 = toi2(locsPos2(idxPos2));
            peaksPosVals2 = pksPos2(idxPos2);
            
            % negative peak T2
            [pksNeg2,locsNeg2,~,promsNeg2] = findpeaks(-vals(toiIdx2,iC));
            [promNeg2, idxNeg2] = sort(promsNeg2,1,'descend');
            idxNeg2 = sort(idxNeg2(1:nPk));
            promNeg2 = promsNeg2(idxNeg2); 
            peaksNeg2 = toi2(locsNeg2(idxNeg2));
            peaksNegVals2 = -pksNeg2(idxNeg2);
         
            % Avg Peak
            peaksPosValsAvg = mean([peaksPosVals1,peaksPosVals2]); 
            peaksNegValsAvg = mean([peaksNegVals1,peaksNegVals2]); 
            promPosAvg = mean([promPos1,promPos2]);
            promNegAvg = mean([promNeg1,promNeg2]); 
            
            % Pos or Neg peak more prominent
            if promPosAvg >  promNegAvg
                promDir = 1; 
                grandPromAvg = promPosAvg; 
                grandPeak1 = peaksPosVals1; % Peak1 val
                grandIdx1 = peaksPos1; % Peak1 idx 
                grandPeak2 = peaksPosVals2; % Peak1 val
                grandIdx2 = peaksPos2; % Peak1 idx 
            elseif promNegAvg >  promPosAvg
                promDir = -1; 
                grandPromAvg = promNegAvg; 
                grandPeak1 = peaksNegVals1; 
                grandIdx1 = peaksNeg1;
                grandPeak2 = peaksNegVals2; 
                grandIdx2 = peaksNeg2;
            end
            
            Pk.peaksPos1(iC,:)       = peaksPos1;
            Pk.peaksPosVals1(iC,:)   = peaksPosVals1;
            Pk.peaksNeg1(iC,:)       = peaksNeg1;
            Pk.peaksNegVals1(iC,:)   = peaksNegVals1;
            
            Pk.peaksPos2(iC,:)       = peaksPos2;
            Pk.peaksPosVals2(iC,:)   = peaksPosVals2;
            Pk.peaksNeg2(iC,:)       = peaksNeg2;
            Pk.peaksNegVals2(iC,:)   = peaksNegVals2;
            
            Pk.promPos1(iC,:)        = promPos1; 
            Pk.promNeg1(iC,:)        = promNeg1; 
            Pk.promPos2(iC,:)        = promPos2; 
            Pk.promNeg2(iC,:)        = promNeg2; 
            
            Pk.peaksPosValsAvg(iC,:) = peaksPosValsAvg; 
            Pk.peaksNegValsAvg(iC,:) = peaksNegValsAvg; 
            Pk.promPosAvg(iC,:)      = promPosAvg; 
            Pk.promNegAvg(iC,:)      = promNegAvg; 
            
            Pk.promDir(iC,:)         = promDir; 
            Pk.grandPromAvg(iC,:)    = grandPromAvg; 
            Pk.grandPeak1(iC,:)      = grandPeak1; 
            Pk.grandIdx1(iC,:)       = grandIdx1; 
            Pk.grandPeak2(iC,:)      = grandPeak2; 
            Pk.grandIdx2(iC,:)       = grandIdx2; 
        end
    end
   
    % sort dirProm
    [sortGrandProm, idxDirProm] = sort(Pk.grandPromAvg,'descend'); 
    % passCh = idxDirProm(1:nChOI); overriding for TANoise 20Hz!!! BUT
    % uncomment otherwise
    
    %% figure: avg ERF per channel, PosNeg peaks T1/T2, toi1/toi2, passCh
    % subplot sizing
    nRow = 7;
    nBreak = 30;
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
        plot(toi1,toiY1,'k','LineWidth',3)
        plot(toi2,toiY2,'Color',[0.5,0.5,0.5],'LineWidth',3)
        plot(t,vals(:,iC),'LineWidth',1,'Color',[0.7 0.7 0.7])
        if sum(ismember(passCh,iC))>0
            if Pk.promDir(iC) == 1
                pcolor = 'b'; 
            elseif Pk.promDir(iC) == -1
                pcolor = 'g';
            end
            plot(t,vals(:,iC),'LineWidth',1,'Color',pcolor)
        end
        plot(Pk.peaksPos1(iC,:), Pk.peaksPosVals1(iC,:), '.b', 'MarkerSize', 10)
        plot(Pk.peaksNeg1(iC,:), Pk.peaksNegVals1(iC,:), '.g', 'MarkerSize', 10)
        plot(Pk.peaksPos2(iC,:), Pk.peaksPosVals2(iC,:), '.b', 'MarkerSize', 10)
        plot(Pk.peaksNeg2(iC,:), Pk.peaksNegVals2(iC,:), '.g', 'MarkerSize', 10)
        
        % prominence lines 
        plot([Pk.peaksPos1(iC,:),Pk.peaksPos1(iC,:)],[Pk.peaksPosVals1(iC,:),Pk.peaksPosVals1(iC,:)-Pk.promPos1(iC,:)],'b')
        plot([Pk.peaksNeg1(iC,:),Pk.peaksNeg1(iC,:)],[Pk.peaksNegVals1(iC,:),Pk.peaksNegVals1(iC,:)+Pk.promNeg1(iC,:)],'g')
        plot([Pk.peaksPos2(iC,:),Pk.peaksPos2(iC,:)],[Pk.peaksPosVals2(iC,:),Pk.peaksPosVals2(iC,:)-Pk.promPos2(iC,:)],'b')
        plot([Pk.peaksNeg2(iC,:),Pk.peaksNeg2(iC,:)],[Pk.peaksNegVals2(iC,:),Pk.peaksNegVals2(iC,:)+Pk.promNeg2(iC,:)],'g')
        
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        yline(0,'--');
        title(sprintf('channel %d, dir %d',iC,Pk.promDir(iC)))
        rd_supertitle2(sprintf('%s: pass channels [%s]', selectionMethod, num2str(passCh'))) 
    end
    
    %% figure: only channels that pass threshold, sorted by max avg prominence 
    figure
    set(gcf,'Position',[100 100 1500 1200])
    hold on
    for iC = 1:nChOI
        subplot(nChOI,1,iC) % one subplot per channel
        iPassCh = passCh(iC); 
        hold on
        if Pk.promDir(iPassCh) == 1
            pcolor = 'b';
        elseif Pk.promDir(iPassCh) == -1
            pcolor = 'g';
        end
        plot(t,vals(:,iPassCh),'LineWidth',1,'Color',pcolor)
        plot(toi1,toiY1,'k','LineWidth',3)
        plot(toi2,toiY2,'Color',[0.5,0.5,0.5],'LineWidth',3)
        plot(Pk.grandIdx1(iPassCh),Pk.grandPeak1(iPassCh),'.k', 'MarkerSize', 10)
        plot(Pk.grandIdx2(iPassCh),Pk.grandPeak2(iPassCh),'.k', 'MarkerSize', 10)
        xlim(xlims)
        vline(p.eventTimes,'k',p.eventNames)
        yline(0,'--');
        title(sprintf('channel %d, dir %d',iPassCh,Pk.promDir(iPassCh)))
        rd_supertitle2(sprintf('%s: pass channels [%s]', selectionMethod, num2str(passCh))) 
    end
    
    %% figure: scatterplot prom and peak 
    figure 
    subplot(2,1,1)
    hold on
    title('T1')
    scatter(Pk.peaksPosVals1,Pk.promPos1,'b')
    scatter(Pk.peaksNegVals1,Pk.promNeg1,'g')
    legend('Pos','Neg')
    xlabel('peak val')
    ylabel('prominence')
    subplot(2,1,2)
    hold on
    title('T2')
    scatter(Pk.peaksPosVals2,Pk.promPos2,'b')
    scatter(Pk.peaksNegVals2,Pk.promNeg2,'g')
    xlabel('peak val')
    ylabel('prominence')
    
    %% figure: prominence and peak overlaid all channels 
    figure
    set(gcf,'Position',[100 100 1200 500])
    subplot(2,1,1)
    hold on
    plot(toi1,toiY1,'k','LineWidth',3)
    plot(toi2,toiY2,'Color',[0.5,0.5,0.5],'LineWidth',3)
    for iC = 1:nCh
        plot(Pk.peaksPos1(iC,:), Pk.peaksPosVals1(iC,:), '.b', 'MarkerSize', 10)
        plot(Pk.peaksNeg1(iC,:), Pk.peaksNegVals1(iC,:), '.g', 'MarkerSize', 10)
        plot(Pk.peaksPos2(iC,:), Pk.peaksPosVals2(iC,:), '.b', 'MarkerSize', 10)
        plot(Pk.peaksNeg2(iC,:), Pk.peaksNegVals2(iC,:), '.g', 'MarkerSize', 10)
        
        % prominence lines
        plot([Pk.peaksPos1(iC,:),Pk.peaksPos1(iC,:)],[Pk.peaksPosVals1(iC,:),Pk.peaksPosVals1(iC,:)-Pk.promPos1(iC,:)],'Color',[0 0 1 0.25])
        plot([Pk.peaksNeg1(iC,:),Pk.peaksNeg1(iC,:)],[Pk.peaksNegVals1(iC,:),Pk.peaksNegVals1(iC,:)+Pk.promNeg1(iC,:)],'Color',[0 1 0 0.25])
        plot([Pk.peaksPos2(iC,:),Pk.peaksPos2(iC,:)],[Pk.peaksPosVals2(iC,:),Pk.peaksPosVals2(iC,:)-Pk.promPos2(iC,:)],'Color',[0 0 1 0.25])
        plot([Pk.peaksNeg2(iC,:),Pk.peaksNeg2(iC,:)],[Pk.peaksNegVals2(iC,:),Pk.peaksNegVals2(iC,:)+Pk.promNeg2(iC,:)],'Color',[0 1 0 0.25])
    end
    xlim(xlims)
    ylabel('amp')
    xlabel('time (ms)')
    vline(p.eventTimes,'k',p.eventNames)
    yline(0,'--');
    title('T1 T2 peaks and prominence')
    
    subplot(2,1,2)
    hold on
    histogram(Pk.grandIdx1,t) 
    histogram(Pk.grandIdx2,t)
    ylabel('count')
    xlabel('time (ms)')
    vline(p.eventTimes,'k',p.eventNames)
    title('grand peak times')
    
    %% figure topo
    figure
    set(gcf,'Position',[100 100 800 800])
    meg_topoplot(Pk.grandPromAvg, [], [], passCh') % avg trial, avg time, highlight channels
    title(sprintf('%s: pass channels [%s]',selectionMethod, num2str(passCh')))
    colorbar
    
    %% store results of peaks analysis
    Pk.t                = t;
    Pk.vals             = vals;
    Pk.nPk              = nPk;
    Pk.sortGrandProm    = sortGrandProm; 
    Pk.idxDirProm       = idxDirProm; 
    Pk.passCh           = passCh;  
    

%% 2. select based on max(abs(amp)) crossing threshhold
case 'maxAmp' 
       
    % plot bar check max(abs(amp)) channels
    figure
    set(gcf,'Position',[100 100 800 300])
    hold on
    b = bar(toiMax,'FaceColor','flat'); 
    b.CData = repmat([0.8 0.8 0.8],[nCh,1]);
    b.CData(passCh,:) = repmat([0 1 0],[nChOI,1]); % highlight channels > threshold 
    hline(thresholdVal,'k',sprintf('threshold = %d',thresholdVal))
    title(sprintf('channels [%s] above threshold',num2str(passCh))) 
    ylabel('max(abs(amp))')
    xlabel('channel')
    
    % plot check group across time
    figure
    hold on 
    plot(t,absVals)
    plot(toi,toiY,'k','LineWidth',3) % highlight black time of interest for checking treshold
    xlim(xlims)
    ylim([0,1.1*max(toiMax)])
    vline(p.eventTimes,'k',p.eventNames)
    hline(thresholdVal,'k',sprintf('threshold = %d',thresholdVal))
    xlabel('time (ms)')
    ylabel('abs(amp)')
    
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
        yline(thresholdVal,'r');
        if ismember(iC,passCh)
            yline(thresholdVal,'g'); % if channel passes threshold, turn line green
        end
        title(sprintf('channel %d',iC))
        rd_supertitle(sprintf('channels [%s] above threshold = %d',num2str(passCh),thresholdVal)) 
    end
   
    % plot only channels that pass threshold, sorted by max 
    figure
    set(gcf,'Position',[100 100 2400 1200])
    hold on
    for iC = 1:nChOI
        subplot(nChOI,1,iC) % one subplot per channel
        hold on
        plot(t,absVals(:,toiMaxChSort(iC)))
        plot(toi,toiY,'k','LineWidth',3)
        xlim(xlims)
        ylim([0,1.1*max(toiMax)])
        vline(p.eventTimes,'k',p.eventNames)
        yline(thresholdVal,'g');
        title(sprintf('channel %s',num2str(toiMaxChSort(iC))))
        rd_supertitle2(sprintf('channels above threshold = %d',thresholdVal)) 
    end

    disp(passCh)
    
    % save threshold info 
    T.passCh = passCh; % channels that pass threshold 
    T.nChOI = nChOI; % number of channels that pass threshold 
    T.thresholdVal = thresholdVal; % threshold value 
    T.thresholdPrctile = thresholdPrctile;  % threshold percentile of channels
    T.toi = toi; % time of interest to look for max(abs(amp))
    T.toiMax = toiMax; % max(abs(amp)) values within toi 
    T.toiMaxSort = toiMaxSort; % sorted toiMax 
    T.toiMaxChSort = toiMaxChSort; % channel idx of sorted toiMax
    T.maxIdx = maxIdx; % idx of max w reference to epoch timing 
    
    %% topo figure
    % set(gcf,'Position',[100 100 800 800])
    % meg_multiplot(vals', [], [], passCh') % ERF all channels avg trial
    figure
    set(gcf,'Position',[100 100 800 800])
    meg_topoplot(toiMax', [], [], passCh') % avg trial, avg time, highlight channels
    title(sprintf('channels [%s] above threshold = %d',num2str(passCh),thresholdVal))
    
    %% avg peaks, treshold
    thresholdVal = prctile(toiMax,thresholdPrctile); % cutoff val
    passCh = find(toiMax > thresholdVal);
    
end

%% save figs

fH = sort(double(findobj(0,'Type','figure')));
figNames = {'ERF_Peak1','ERF_Peak2','ERF_Peak3','ERF_Peak4','ERF_Peak5','ERF_Peak6',...
            'ERF_passCh','scatterPromPeak','histogramPeakTiming','topo_avgProm'};


% if selectThreshold
%     figDirThresh = sprintf('%sthreshold%d',figDir,nChOI); 
%     figNames = {'histMaxAmp','ERF_ThreshAllChannels','ERF_Thresh1','ERF_Thresh2','ERF_Thresh3','ERF_Thresh4','ERF_Thresh5','ERF_passCh','topoPassCh'};
%     if ~exist(figDirThresh,'dir')
%         mkdir(figDirThresh)
%     end
%     rd_saveAllFigs(fH, figNames, sessionDir, figDirThresh, [])
% end

% figNames = {'TopChannels_Alpha'};
% figNames = {'ERFPeaks1','ERFPeaks2','ERFPeaks3','ERFPeaks4','ERFPeaks5','TopoERF'}
rd_saveAllFigs(fH, figNames, sessionDir, figDir, [])

%% save info 

save(sprintf('%s/Pk.mat',matDir),'Pk')
% save(sprintf('%s/C.mat',matDir),'C')

% if selectThreshold
%     save(sprintf('%s/T.mat',matDir),'T') % top channels by thresholding
% end


end



