function [A,fH,figNames] = meg_permutationTest(data1,data2,nPerm,thresh,stat)

% MEG_PLOTERP(data,selectedChannels,plotSingleTrial,plotAvgTrial)
% 
% implements nonparametric statistical testing, see Maris & Oostenveld, 2007: https://www.researchgate.net/publication/6316066_Nonparametric_statistical_testing_of_EEG-_and_MEG-data
%   - ttest data1, data2
%   - randomly permutes data into 2 conds nPerm times 
%   - ttest on random permutations
%   - rd_clusterStat: uses threshold to determine max cluster of tvals
%   - compares distribution of max cluster stats from permutations to true
%   data 
%
% INPUT
% data1, data2
%   data should be 3 dimensional (e.g. freq x time x subject) 
%   data by cond (between 2 conditions within subjects) 
%   or data by group (bewteen 2 subject groups) 
% nPerm
%   number of random permutations 
% thresh 
%   statistical threshold, eg ts >.8;
% stat 
%   cluster stat method, either 'sum' or 'mean'

% OUPUT
% A
%   s
% fH
%   figure handle
% figNames 
%   figure names 
%
% Karen Tian
% November 2020

%% 

hAll = []; pvalAll = []; 
clusterStats = []; maxAbsClusterStatAll = []; 

%% combine data

sizeData1 = size(data1); 
sizeData2 = size(data2); 
dimension = size(sizeData1,2); 

dataAll = cat(dimension,data1,data2); 

%% ttest on true data (data1, data2)

h_true = zeros(size(im,1),size(im,2));
pval_true = zeros(size(im,1),size(im,2));
tval_true = zeros(size(im,1),size(im,2));
for r = 1:sizeData1(1)
    for c = 1:sizeData1(2)
        if all(~isnan(data1(r,c,:))) && all(~isnan(data2(r,c,:)))
            [h_true(r,c),pval_true(r,c),~,stats] = ttest2(data1(r,c,:),data2(r,c,:));
            tval_true(r,c) = stats.tstat;
        else
            h_true(r,c) = NaN;
            pval_true(r,c) = NaN;
        end
    end
end

% thresh = im>0; 
% [clusterStats, maxAbsClusterStat, C] = rd_clusterStat(imabs,thresh);

%% 

% initialize variables 
nConds = size(dataAll,3); 
clusterStats = cell(1,nPerm); 
maxAbsClusterStat = NaN(1,nPerm); 
hAll = []; 

% check number of conds to randomly permute, if even split, if uneven
% randomly assign last odd 
if rem(nConds,2)==0 % even
    numConds = 'isEven'; 
else
    numConds = 'isOdd'; 
end

for i = 1:nPerm % number of random shuffles 
    idxPerm = randperm(nConds);
    switch numConds
        case 'isEven'
            idx1 = idxPerm(1:nConds/2);
            idx2 = idxPerm(nConds/2+1:nConds); 
        case 'isOdd'
            idx1 = idxPerm(1:(nConds-1)/2);
            idx2 = idxPerm((nConds+1)/2:nConds);
    end
    
    data1_rand = dataAll(:,:,idx1); 
    data2_rand = dataAll(:,:,idx2);
    
    % point by point ttest (or ranksum for Wilcoxon rank test? equivalent to
    % Mann-Whitney U-test for unequal samples)
    
    tic
    h_rand = NaN(size(im,1),size(im,2));
    pval_rand = NaN(size(im,1),size(im,2));
    tval_rand = NaN(size(im,1),size(im,2));
    for r = 1:sizeData1(1)
        for c = 1:sizeData1(2)
            if all(~isnan(data1_rand(r,c,:))) && all(~isnan(data2_rand(r,c,:)))
                [h_rand(r,c),pval_rand(r,c),~,stats] = ttest2(data1_rand(r,c,:),data2_rand(r,c,:));
                tval_rand(r,c) = stats.tstat;
%             else
%                 h_rand(r,c) = NaN;
%                 pval_rand(r,c) = NaN;
            end
        end
    end
    toc
    
    tic
    thresh = abs(tval_rand)>2.3646; % edit w appropriate critical tval
    [clusterStats{i}, maxAbsClusterStat(i),~] = rd_clusterStat(tval_rand, thresh, stat);
    toc

    % save iteration statistics 
    hAll = cat(3,hAll,h_rand);
    % save tval_rand all 
    % pvalAll = cat(3,pvalAll,pval);
    % clusterStatsAll = cat(1,clusterStatsAll,clusterStats); 
    disp(i) % stopped 522 
    disp(maxAbsClusterStat(i)) 
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

titleText = 'permutation_maxAbsClusterStat'; 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))


%% plot 1st 50 permuted hypothesis and pvals 

figure
set(gcf, 'Position',  [100, 100, 2000, 800])
nfig = 1; 
xtick = [2001,2300+2001]; 
for i = 1:50
    hold on
    subplot(5,10,i)
    imagesc(hAll(:,:,i))
    xlim(xlims)
    ylim(ylims)
    % colorbar
    % caxis([0,1])
    
    set(gca,'YDir','normal')
    set(gca,'XTick',xtick)
    set(gca,'YTick',ytick)
    set(gca,'XTickLabel',toi(xtick))
    set(gca,'YTickLabel',foi(ytick))
    
    hold on
    for iEv = 1:numel(p.eventTimes)
        vline(p.eventTimes(iEv)+abs(p.tstart),'w');
    end
    % meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    %     title(sprintf('%s ITPC',subjectNames{i}))
    %     titleText = sprintf('ITPC_%s',subjectNames{i});
    % print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
    nfig = nfig + 1;
    if nfig > 10
        nfig = 1;
    end
    if i == 1
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')

    end
end

rd_supertitle2('ITPC permutation hypothesis tests')
titleText = 'ITPCpermutation_clusterStat_histogram_normalized'; % 'ITPCpermutation'; 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% true data difference

data1group = nanmean(data1,3);
data2group = nanmean(data2,3);

im = data1group - data2group; 

tval = NaN(size(im,1),size(im,2));
h = NaN(size(im,1),size(im,2));
pval = NaN(size(im,1),size(im,2));
for r = 1:size(im,1)
    for c = 1:size(im,2)
        if all(~isnan(data1(r,c,:))) && all(~isnan(data2(r,c,:)))
            [h(r,c),pval(r,c),~,stats] = ttest2(data1(r,c,:),data2(r,c,:));
            tval(r,c) = stats.tstat;
        end
    end
end

clusterStats_true = []; maxAbsClusterStat_true = []; 
thresh = abs(tval)>2.3646; % edit w appropriate critical tval
[clusterStats_true, maxAbsClusterStat_true,C] = rd_clusterStat(tval, thresh, stat);
sigIm = zeros(size(im)); 
for i = 1:numel(C.PixelIdxList)
    sigIm(C.PixelIdxList{i}) = 1;
    if clusterStats_true(i) > permuteThresh 
        sigIm(C.PixelIdxList{i}) = 10; 
    end
end

%% difference figure hypothesis 

figure
set(gcf, 'Position',  [100, 100, 1000, 600])
xtick = [2001,2300+2001]; 

hold on
imagesc(h)
xlim(xlims)
ylim(ylims)
% colorbar
% caxis([0,1])

set(gca,'YDir','normal')
set(gca,'XTick',xtick)
set(gca,'YTick',ytick)
set(gca,'XTickLabel',toi(xtick))
set(gca,'YTickLabel',foi(ytick))

hold on
for iEv = 1:numel(p.eventTimes)
    vline(p.eventTimes(iEv)+abs(p.tstart),'w');
end
% meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
%     title(sprintf('%s ITPC',subjectNames{i}))
%     titleText = sprintf('ITPC_%s',subjectNames{i});
% print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
xlabel('Time (s)')
ylabel('Frequency (Hz)')

rd_supertitle2('hypothesis test true data')
titleText = 'ITPCspectrogram_trueData_hypothesistest'; % 'ITPCpermutation';
% print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% plot histogram of max abs permutation test and find 95 percentile cutoff + true data histogram
figure
set(gcf, 'Position',  [100, 100, 800, 300])
hold on 
histogram(maxAbsClusterStat,nPerm/20) 
histogram(abs(clusterStats_true)) 
permuteThresh = prctile(maxAbsClusterStat,95); 
permuteThresh90 = prctile(maxAbsClusterStat,90); 
xline(permuteThresh,'k',permuteThresh,'LineWidth',2)
xline(permuteThresh90,'k',permuteThresh90,'LineWidth',2)
xline(maxAbsClusterStat_true,'r',maxAbsClusterStat_true,'LineWidth',2)

set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off
xlabel('Max abs cluster based permuatation sum statistic')
ylabel('Count')

titleText = 'permutation_maxAbsClusterStat_realData_histogram'; 
% print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))


