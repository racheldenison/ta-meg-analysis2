function [A] = meg_nonparametricTest(data1,data2,nPerm,stat,variableName,condName,t)
% MEG_NONPARAMETRICTEST()
% 
% performs nonparametric testing on time series data with two conditions (Maris & Oostenveld 2007)
%   - ttest data1, data2
%   - randomly permutes data into 2 conds nPerm times 
%   - ttest on random permutations
%   - rd_clusterStat: uses threshold to determine max cluster of tvals
%   - compares distribution of max cluster stats from permutations to true
%   data 
% 
%
% INPUT
% data1, data2
%   data should be 3 dimensional (e.g. freq x time x subject) 
%       or 2 dimensional (e.g. time x subject) 
%   data by cond (between 2 conditions within subjects) 
%   or data by group (bewteen 2 subject groups) 
% nPerm
%   number of random permutations 
% thresh 
%   statistical threshold, eg ts >.8; 
% stat 
%   cluster stat method, either 'sum' or 'mean'
% variable name 
%   saves figures w appropriate variable name 
% t 
%   time vector 

% OUPUT
% A
%   shuffled data
% fH
%   figure handle
% figNames 
%   figure names 
%
% Karen Tian
% April 2021 
% Updated Oct 2022 

%% check inputs 
% for ITPC permutation test, load groupA_ITPCspectrogram_byAtt.mat
if nargin<7 
    t = size(data1,1); 
end
if nargin<6
    condName = ''; 
end
if nargin<5
    variableName = ''; 
end
if nargin<4
    stat = 'sum';
    disp('statistic not specific, default sum')
end
if nargin<3 
    nPerm = 1000; 
    disp('number of permutations not specified, default 1000')
end
if nargin<2
    error('data matrix of 2nd condition not specified') 
end
if nargin<1
    error('data matrix of 1st condition not specified') 
end

%% initialize 
hAll = []; pvalAll = []; 
clusterStats = []; maxAbsClusterStatAll = []; 

plotFig = 1; 
saveFig = 1; 
figDir = sprintf('%s/permutationTest',pwd); 
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% combine data
sizeData1 = size(data1); 
sizeData2 = size(data2); 
dimension = size(sizeData1,2); 

dataAll = cat(dimension,data1,data2); % combine data from cond 1 and 2

df = sizeData1(dimension)-1; % -1 for 1 sample ttest 
tail = 'two'; % 'one' 'two'
criticalTVal = CritT(0.05,df,tail); 

%% parameters 
p = meg_params('TANoise_ITPCsession8'); 
shuffledVal = []; 

% t = p.tstart:p.tstop; 

expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

%% data
data = cat(dimension+1,data1,data2); 

% vals_cueT1 = data(:,:,1); 
% vals_cueT2 = data(:,:,2); 

nSamples = size(data,1);
xlims = [1,sizeData1(dimension-1)]; 

random_data = repmat(data,[1 1 1 nPerm]);  % t x subject x cond x perm

%% permute, shuffle attention condition 
% figure
% set(gcf, 'Position',  [100, 100, 2000, 200])
% if nPerm > 10
%     nPlot = 5;
% else
%     nPlot = nPerm;
% end
% for i = 1:nPerm
%     % permute cue conditions for these indices
%     idxRand1 = logical(randi([0 1],1,size(data,2)));
% %     idxRand2 = idxRand1+1;
% %     idxRand2(idxRand2==2) = 0;
% %     idxRand2 = logical(idxRand2);
%     
%     random_data(:,idxRand1,1,i) = data(:,idxRand1,2);
%     random_data(:,idxRand1,2,i) = data(:,~idxRand1,1);
%     
%     if plotFig && i<nPlot 
%         subplot (1,nPlot+1,1)
%         hold on
%         plot(data(:,:,1),'Color',p.cueColors(1,:))
%         plot(data(:,:,2),'Color',p.cueColors(2,:))
%         title('true data')
%         xlim([1 numel(t)])
%         set(gca,'TickDir','out');
%         ax = gca;
%         ax.LineWidth = 1.5;
%         ax.XColor = 'black';
%         ax.YColor = 'black';
%         ax.FontSize = 12;
%         box off
%         
%         subplot (1,nPlot+1,i+1)
%         hold on
%         plot(random_data(:,:,1,i),'Color',p.cueColors(1,:))
%         plot(random_data(:,:,2,i),'Color',p.cueColors(2,:))
%         title('random data')
%         xlim([1 numel(t)])
%         set(gca,'TickDir','out');
%         ax = gca;
%         ax.LineWidth = 1.5;
%         ax.XColor = 'black';
%         ax.YColor = 'black';
%         ax.FontSize = 12;
%         box off
%     end 
% end
% titleText = sprintf('20Hz_%s_randAttCond_%d,%d',variableName,t(1),t(end));
% if saveFig 
%     print(gcf,'-dpng','-painters',titleText)
% end

%% plot difference
% if plotFig
%     if nPerm > 10 
%         nPlots = 10; 
%     else
%         nPlots = nPerm; 
%     end
%     figure
%     set(gcf, 'Position',  [100, 100, 2000, 200])
%     for i = 1:nPlots
%         subplot (1,nPlots+1,1)
%         hold on
%         plot(data(:,:,1)-data(:,:,2),'Color',p.cueColors(3,:))
%         title('true data')
%         xlim([1 numel(t)])
%         set(gca,'TickDir','out');
%         ax = gca;
%         ax.LineWidth = 1.5;
%         ax.XColor = 'black';
%         ax.YColor = 'black';
%         ax.FontSize = 12;
%         box off
%         
%         subplot (1,nPlots+1,i+1)
%         hold on
%         plot(random_data(:,:,1,i)-random_data(:,:,2,i),'Color',p.cueColors(3,:))
%         title('random data')
%         xlim([1 numel(t)])
%         set(gca,'TickDir','out');
%         ax = gca;
%         ax.LineWidth = 1.5;
%         ax.XColor = 'black';
%         ax.YColor = 'black';
%         ax.FontSize = 12;
%         box off
%     end
% end
% titleText = sprintf('20Hz_%s_randCond_cueT1-cueT2_%d,%d',variableName,t(1),t(end));
% if saveFig
%     print(gcf,'-dpng','-painters',titleText)
% end

%% plot data permutation 
% figure
% set(gcf, 'Position',  [100, 100, 800, 200])
% for i = 1:20
%     hold on
%     pLine = plot(shuffledVal(:,i));  
%     pLine.Color = [p.cueColors(3,:) 0.4]; 
% end
% plot(data,'LineWidth',2,'Color',p.cueColors(3,:))
% xlim([0 105])
% 
% set(gca,'TickDir','out');
% ax = gca;
% ax.LineWidth = 1.5;
% ax.XColor = 'black';
% ax.YColor = 'black';
% ax.FontSize = 12;
% box off
% xlabel('Time')
% ylabel('Amplitude (T)')
% title('permuted time series (nPerm = 1000)')
% 
% titleText = 'PermutedTS_cueT1-cueT2';
% if saveFig
%     print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
% end

%% ttest on true data (data1, data2) uncorrected 
if dimension == 3
    h_true = NaN(sizeData1(1),sizeData1(2));
    pval_true = NaN(sizeData1(1),sizeData1(2));
    tval_true = NaN(sizeData1(1),sizeData1(2));
    for r = 1:sizeData1(1)
        for c = 1:sizeData1(2)
            if all(~isnan(data1(r,c,:))) && all(~isnan(data2(r,c,:)))
                [h_true(r,c),pval_true(r,c),~,stats] = ttest(squeeze(data1(r,c,:)-data2(r,c,:)));
                tval_true(r,c) = stats.tstat;
            else
                h_true(r,c) = NaN;
                pval_true(r,c) = NaN;
            end
        end
    end
    
elseif dimension == 2
    h_true = NaN(sizeData1(1),1);
    pval_true = NaN(sizeData1(1),1);
    tval_true = NaN(sizeData1(1),1);
    for r = 1:sizeData1(1) % time
        if all(~isnan(data1(r,:))) && all(~isnan(data2(r,:)))
            [h_true(r),pval_true(r),~,stats] = ttest(data1(r,:),data2(r,:));
            tval_true(r) = stats.tstat;
        else
            h_true(r) = NaN;
            pval_true(r) = NaN;
        end
    end
end

thresh = abs(tval_true)>=criticalTVal; 
[clusterStats_true, maxAbsClusterStat_true,C] = rd_clusterStat(tval_true, thresh, stat);

% PLOT 
if dimension==3
    figure
        imagesc(h_true)  
    figure
        imagesc(mean(data1,dimension,'omitnan')-mean(data2,dimension,'omitnan'))
elseif dimension<3
    figure   
    subplot 211
    hold on
    imagesc(h_true')
    % colorbar
    caxis([0 1])
    xlim(size(h_true'))
    xticks([])
    yticks([])
    meg_figureStyle
    title('point wise ttest (uncorrected), yellow > criticalT')
    
    subplot 212
    hold on
    plot(p.t(t),mean(data1,2,'omitnan'))
    plot(p.t(t),mean(data2,2,'omitnan'))
    meg_figureStyle
    xlim([p.t(t(1)) p.t(t(end))])
    xlabel('Time (ms)')
    ylabel(sprintf('%s',und2space(variableName)))
end
sgtitle(sprintf('%s\n%s\nTrue max abs cluster stat = %0.2f',und2space(variableName),und2space(condName),maxAbsClusterStat_true))

titleText = sprintf('uncorrected_ttest_%s_%s_n%d',variableName,condName,size(data,2));
if saveFig
    print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
end

%% randomly permute conds 
% THIS WAS FOR TESTING UP DONWERS
% if dimension==3
%     nConds = size(dataAll,3); 
% elseif dimension==2
%     % nConds = size(dataAll,2); 
%     nConds = 2; 
% end
% clusterStats = cell(1,nPerm); 
% maxAbsClusterStat = NaN(1,nPerm); 
% hAll = []; 
% tValRandAll = []; 
% 
% % check number of conds to randomly permute
% % if even split, if uneven randomly assign last odd 
% if rem(nConds,2)==0 % even
%     numConds = 'isEven'; 
% else
%     numConds = 'isOdd'; 
% end
% 
% for i = 1:nPerm 
%     idxPerm = randperm(nConds);
%     switch numConds
%         case 'isEven'
%             idx1 = idxPerm(1:nConds/2);
%             idx2 = idxPerm(nConds/2+1:nConds); 
%         case 'isOdd'
%             idx1 = idxPerm(1:(nConds-1)/2);
%             idx2 = idxPerm((nConds+1)/2:nConds);
%     end
%     
%     if dimension==3
%         data1_rand = dataAll(:,:,idx1);
%         data2_rand = dataAll(:,:,idx2);
%     elseif dimension==2
%         data1_rand = dataAll(:,idx1);
%         data2_rand = dataAll(:,idx2);
%     end
%     
%     h_rand=[]; pval_rand=[]; tval_rand=[];
%     % point by point ttest on rand permutation 
%     if dimension==3
%         tic
%         h_rand = NaN(sizeData1(1),sizeData1(2));
%         pval_rand = NaN(sizeData1(1),sizeData1(2));
%         tval_rand = NaN(sizeData1(1),sizeData1(2));
%         for r = 1:sizeData1(1)
%             for c = 1:sizeData1(2)
%                 if all(~isnan(data1_rand(r,c,:))) && all(~isnan(data2_rand(r,c,:)))
%                     [h_rand(r,c),pval_rand(r,c),~,stats] = ttest(squeeze(data1_rand(r,c,:)-data2_rand(r,c,:)));
%                     tval_rand(r,c) = stats.tstat;
%                 end
%             end
%         end
%         toc
%     elseif dimension==2
%         tic
%         h_rand = NaN(sizeData1(1),1); % just rows
%         pval_rand = NaN(sizeData1(1),1);
%         tval_rand = NaN(sizeData1(1),1);
%         for r = 1:sizeData1(1)   
%             if all(~isnan(data1_rand(r,:))) && all(~isnan(data2_rand(r,:)))
%                 [h_rand(r),pval_rand(r),~,stats] = ttest(squeeze(data1_rand(r,:)-data2_rand(r,:)));
%                 tval_rand(r) = stats.tstat;
%             end
%         end
%         toc
%     end
%     
%     tic
%     thresh = abs(tval_rand)>criticalTVal; % edit w appropriate critical tval
%     [clusterStats{i}, maxAbsClusterStat(i),~] = rd_clusterStat(tval_rand, thresh, stat);
%     toc
% 
%     % save iteration statistics 
%     hAll = cat(3,hAll,h_rand);
%     tValRandAll = cat(3,tValRandAll,tval_rand); 
%     % save tval_rand all 
%     % pvalAll = cat(3,pvalAll,pval);
%     % clusterStatsAll = cat(1,clusterStatsAll,clusterStats); 
%     disp(i) % stopped 522 
%     disp(maxAbsClusterStat(i)) 
% end

%% randomly permute conds 
% NEW SHUFFLE CUE COND AT SUBJECT LEVEL 

% if plotFig
%     if nPerm < 10
%         nPlots = nPerm;
%     else
%         nPlots = 10;
%     end
%     count = 0;
%     figure
%     set(gcf, 'Position',  [100, 100, 600, 1200])
% end

clusterStats = cell(1,nPerm); 
maxAbsClusterStat = NaN(1,nPerm); 
hAll = []; 
tValRandAll = []; 

for i = 1:nPerm 
    idxRand = logical(randi([0 1],1,size(data,2))); % if 1 shuffle, if 0 do not shuffle 
    dataRand = NaN(size(data)); % time x subject x cond 
    for iS = 1:size(data1,2)
        if idxRand(iS) % swap 
            dataRand(:,iS,2) = data(:,iS,1); 
            dataRand(:,iS,1) = data(:,iS,2); 
        elseif ~idxRand(iS) % don't swap 
            dataRand(:,iS,1) = data(:,iS,1); 
            dataRand(:,iS,2) = data(:,iS,2); 
        end
    end
    
    data1_rand = dataRand(:,:,1);
    data2_rand = dataRand(:,:,2);

    h_rand=[]; pval_rand=[]; tval_rand=[];
    % point by point ttest on rand permutation 
    tic
    h_rand = NaN(sizeData1(1),1); % just rows
    pval_rand = NaN(sizeData1(1),1);
    tval_rand = NaN(sizeData1(1),1);
    for r = 1:sizeData1(1)
        if all(~isnan(data1_rand(r,:))) && all(~isnan(data2_rand(r,:)))
            [h_rand(r),pval_rand(r),~,stats] = ttest(data1_rand(r,:),data2_rand(r,:));
            tval_rand(r) = stats.tstat;
        end
    end
    toc
    
    tic
    thresh = abs(tval_rand)>criticalTVal; % edit w appropriate critical tval
    [clusterStats{i}, maxAbsClusterStat(i),~] = rd_clusterStat(tval_rand, thresh, stat);
    toc

    % save iteration statistics 
    hAll = cat(3,hAll,h_rand);
    tValRandAll = cat(3,tValRandAll,tval_rand); 
    % save tval_rand all 
    % pvalAll = cat(3,pvalAll,pval);
    % clusterStatsAll = cat(1,clusterStatsAll,clusterStats); 
    disp(i) % stopped 522 
    disp(maxAbsClusterStat(i)) 
    
%     if plotFig
%         if i<=nPlots
%             subplot (nPlots,2,i*2-1)
%             hold on
%             imagesc(h_rand')
%             xlim(size(h_rand'))
%             xticks([])
%             yticks([])
%             meg_figureStyle
%             
%             subplot (nPlots,2,i*2)
%             hold on
%             plot(p.t(t),mean(data1_rand,2,'omitnan'))
%             plot(p.t(t),mean(data2_rand,2,'omitnan'))
%             meg_figureStyle
%             xlim([p.t(t(1)),p.t(t(end))])
%         end
%     end
%     sgtitle(sprintf('%s\n%s\nTrue max abs cluster stat = %0.2f',und2space(variableName),und2space(condName),maxAbsClusterStat_true))
%     titleText = sprintf('perm%d_ttest_%s_%s_n%d',nPlots,variableName,condName,size(data,2));
%     if saveFig
%         print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
%     end
end
    

%% plot histogram of max abs permutation test and find 95 percentile cutoff 
figure
set(gcf, 'Position',  [100, 100, 800, 300])
histogram(maxAbsClusterStat,nPerm/20) 
permuteThresh = prctile(maxAbsClusterStat,95); 
xline(permuteThresh,'k',permuteThresh,'LineWidth',2)

set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off
xlabel('Max abs cluster based permuatation sum statistic')
ylabel('Count')

titleText = sprintf('%s_permutation_maxAbsClusterStat',variableName); 
if saveFig
    print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
end

%% plot permuted vals and hypothesis 


%% plot 1st 50 permuted hypothesis and pvals 
% if dimension==3
%     figure
%     set(gcf, 'Position',  [100, 100, 2000, 800])
%     nfig = 1;
%     xtick = [1,p.eventTimes(3)-p.eventTimes(2)+1];
%     for i = 1:20
%         hold on
%         subplot(2,10,i)
%         imagesc(hAll(:,:,i))
%         xlim(xlims)
%         
%         set(gca,'YDir','normal')
%         set(gca,'XTick',xtick)
%         set(gca,'YTick',ytick)
%         set(gca,'XTickLabel',toi(xtick))
%         set(gca,'YTickLabel',foi(ytick))
%         
%         nfig = nfig + 1;
%         if nfig > 10
%             nfig = 1;
%         end
%         if i == 1
%             xlabel('Time (s)')
%             ylabel('Frequency (Hz)')
%             
%         end
%     end
% elseif dimension==2
%     figure
%     set(gcf, 'Position',  [100, 100, 2000, 500])
%     nfig = 1;
%     xtick = [1,p.eventTimes(3)-p.eventTimes(2)+1];
%     for i = 1:20
%         hold on
%         subplot(2,10,i)
%         plot(t,squeeze(hAll(:,:,i)))
%         xlim(xlims)
%         meg_figureStyle
%         
%         nfig = nfig + 1;
%         if nfig > 10
%             nfig = 1;
%         end
%         if i == 1
%             xlabel('Time (s)')
%             ylabel(sprintf('%s',variableName))
%         end
%     end
% end
% 
% rd_supertitle2(sprintf('%s permutation hypothesis tests',variableName))
% titleText = sprintf('%s_permutation_clusterStat_histogram',variableName); % 'ITPCpermutation'; 
% if saveFig 
%     print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
% end

%% true data difference
if dimension==3 % time x freq 
    data1group = nanmean(data1,3);
    data2group = nanmean(data2,3);
    im = data1group - data2group;
    sigIm = NaN(size(im));
    for i = 1:numel(C.PixelIdxList)
        if abs(clusterStats_true(i)) > permuteThresh
            sigIm(C.PixelIdxList{i}) = .8;
        end
    end
    
    figure
    set(gcf, 'Position',  [100, 100, 1000, 600])
    
    hold on
    imagesc(im)
    xlim(xlims)
    % ylim(ylims)
    % colorbar
    % caxis([0,1])
    
    set(gca,'YDir','normal')
    set(gca,'XTick',xtick)
    set(gca,'YTick',ytick)
    set(gca,'XTickLabel',toi(xtick))
    set(gca,'YTickLabel',foi(ytick))
    
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
elseif dimension==2 %  time (att) 
    data1group = nanmean(data1,2); % average subjects
    data2group = nanmean(data2,2);
    im = data1group - data2group;
    sigIm = NaN(size(im));
    for i = 1:numel(C.PixelIdxList)
        if abs(clusterStats_true(i)) > permuteThresh
            sigIm(C.PixelIdxList{i}) = .8;
        end
    end
    figure
    set(gcf, 'Position',  [100, 100, 1000, 600])
    hold on
    plot(im)
    xlim(xlims)
    xlabel('Time (ms)')
    meg_figureStyle
end

rd_supertitle2('true data difference')
titleText = 'trueData_difference'; % 'ITPCpermutation';
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% true data difference hypothesis 
if dimension==3
    data1group = nanmean(data1,3);
    data2group = nanmean(data2,3);
    im = data1group - data2group;
    
    figure
    set(gcf, 'Position',  [100, 100, 1000, 600])
    
    hold on
    imagesc(im)
    
    redIm = flipud(cat(3,ones(size(sigIm)),zeros(size(sigIm)),zeros(size(sigIm))));
    redImHandle = imagesc(redIm);
    set(redImHandle,'AlphaData',h_true)
    
    xlim(xlims)
    ylim(ylims)
    colorbar
    % caxis([0,1])
    
    set(gca,'YDir','normal')
    set(gca,'XTick',xtick)
    set(gca,'YTick',ytick)
    set(gca,'XTickLabel',toi(xtick))
    set(gca,'YTickLabel',foi(ytick))
    
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    
    rd_supertitle2('true data difference hypothesis')
    titleText = 'trueData_difference_hypothesis'; % 'ITPCpermutation';
    print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
    
    %% difference figure hypothesis
    figure
    set(gcf, 'Position',  [100, 100, 1000, 600])
    
    hold on
    imagesc(im)
    xlim(xlims)
    ylim(ylims)
    colorbar
    % caxis([0,1])
    
    hold on
    redIm = cat(3,ones(size(sigIm)),zeros(size(sigIm)),zeros(size(sigIm)));
    redImHandle = imagesc(redIm);
    set(redImHandle,'AlphaData',sigIm)
    
    set(gca,'YDir','normal')
    set(gca,'XTick',xtick)
    set(gca,'YTick',ytick)
    set(gca,'XTickLabel',toi(xtick))
    set(gca,'YTickLabel',foi(ytick))
    
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    
    rd_supertitle2('hypothesis test true data')
    titleText = 'trueData_hypothesistest'; % 'ITPCpermutation';
    print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
end

%% plot histogram of max abs permutation test and find 95 percentile cutoff + true data histogram
figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on 
histRand = histogram(maxAbsClusterStat,nPerm/20);
histRand.BinWidth = 100; 
histTrue = histogram(abs(clusterStats_true)); 
histTrue.BinWidth = 100; 
permuteThresh = prctile(maxAbsClusterStat,95); 
permuteThresh90 = prctile(maxAbsClusterStat,90); 
xline(permuteThresh,'k',sprintf('95 percentile permutation threshold = %0.2f',permuteThresh),'LineWidth',2)
xline(permuteThresh90,'k',sprintf('90 percentile permutation threshold = %0.2f',permuteThresh90),'LineWidth',2)
xline(maxAbsClusterStat_true,'r',sprintf('true max abs cluster = %0.2f',maxAbsClusterStat_true),'LineWidth',2)

set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off
xlabel('Max abs cluster based permuatation sum statistic')
ylabel('Count')

titleText = sprintf('%s_permutation_maxAbsClusterStat_realData_histogram',variableName); 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% save analysis 
A.cond1 = 'PrecueT1'; % ITPC uppers 
A.cond2 = 'PrecueT2'; % ITPC downers 
A.toi = xlims(1)+p.eventTimes(2):xlims(2)+p.eventTimes(2); 

A.df = df; 
A.criticalTVal = criticalTVal; 
A.h_true = h_true; 
A.pval_true = pval_true; 
A.tval_true = tval_true; 

A.clusterStats_true = clusterStats_true; 
A.maxAbsClusterStat_true = maxAbsClusterStat_true; 
A.maxAbsClusterStat_percentile = invprctile(maxAbsClusterStat,maxAbsClusterStat_true); 
A.C = C; 

A.diffIm = im; 
A.sigDiffIm = sigIm; 

A.hAll_rand = hAll; 
A.tValRandAll = tValRandAll;
A.clusterStats_rand = clusterStats; 
A.maxAbsClusterStatAll_rand = maxAbsClusterStat; 
A.permuteThresh95_rand = permuteThresh; 

save(sprintf('%s_spectrogram_permutationTest_T1_T2_average.mat',variableName),'A','-v7.3')

%% 


