% November 3, 2020 
% kt_session_makeSubject_upDown_group.m
% load 

%% setup 
A.variable = 'power_avgTrial_cueT2'; % 'ITPC_singleTrial_cueT2'; 
A.normalization_type = 'baselineSubtraction_500-1000'; % 'baselineSubtraction_-500-0'; 'none'
A.session = vals_cueT2; % valsAll% groupITPCspectrogram; 

p = meg_params('TANoise_ITPCsession8');
t = p.tstart:p.tstop;
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = p.tstart:1:p.tstop;
tfAmps = [];

%% session --> subject 
count = 1; 
for i = 1:2:19
    vals1 = []; vals2 = [];
    vals1 = A.session(:,:,i); 
    vals2 = A.session(:,:,i+1); 
    vals12 = cat(3,vals1,vals2); 
    meanVals12 = nanmean(vals12,3); 
    A.subject(:,:,count) = meanVals12; 
    count = count+1; 
end

%%  subject --> up/down 
expt = 'TANoise';
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 
up = []; down = []; 
for i = 1:numel(subjectNames)
    vals = []; 
    if ITPCsubject(i) == 1
        vals = A.subject(:,:,i); 
        up = cat(3,up,vals); 
    elseif ITPCsubject(i) == -1 
        vals = A.subject(:,:,i); 
        down = cat(3,down,vals);   
    end
end
meanUp = nanmean(up,3); 
meanDown = nanmean(down,3); 
A.up = up; 
A.down = down; 
A.upMean = meanUp; 
A.downMean = meanDown; 

%% subject --> group
group = []; 
for i = 1:numel(subjectNames)
    vals = []; 
    if ITPCsubject(i) ~= 0
        vals = A.subject(:,:,i); 
        group = cat(3,group,vals); 
    end
end
group = nanmean(group,3); 
A.group = group; 

%% baseline normalization ITPC 
% 500-1000ms 

idxBaseline = abs(p.tstart)+500:abs(p.tstart)+1000; 
A.normMethod = 'baseline_mean_subtraction'; 
A.normBaselineIdx = idxBaseline; 
A.normToi = 500:1000; 

normSubject = []; 
for i = 1:numel(subjectNames)
    vals = []; 
    vals = A.subject(:,idxBaseline,i); 
    meanBaseline = nanmean(vals,2); 
    normVals = A.subject(:,:,i) - meanBaseline; 
    normSubject = cat(3,normSubject,normVals); 
end
A.normSubject = normSubject; 


%% normalized up down group 
up = []; down = []; 
for i = 1:numel(subjectNames)
    vals = []; 
    if ITPCsubject(i) == 1
        vals = A.normSubject(:,:,i); 
        up = cat(3,up,vals); 
    elseif ITPCsubject(i) == -1 
        vals = A.normSubject(:,:,i); 
        down = cat(3,down,vals);   
    end
end
meanUp = nanmean(up,3); 
meanDown = nanmean(down,3); 
A.normUp = meanUp; 
A.normDown = meanDown; 

%% plot normalized subject 
ylims = [1,50]; 
xlims = [size(toi,1),size(toi,2)]; 
toi = p.tstart:1:p.tstop;
ytick = 10:10:numel(foi);
xtick = 1:1000:numel(toi);

for i = 1 % :numel(subjectNames)
    figure
    set(gcf, 'Position',  [100, 100, 400, 300])
    hold on
    imagesc(A.normSubject(:,:,i))
    xlim(xlims)
    ylim(ylims)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar
    caxis([-0.3,0.3])

    set(gca,'YDir','normal')
     set(gca,'XTick',xtick)
     set(gca,'XTickLabel',toi(xtick))
     set(gca,'YTick',ytick)
     set(gca,'YTickLabel',foi(ytick))
     hold on
    for iEv = 1:numel(p.eventTimes)
        vline(p.eventTimes(iEv)+abs(p.tstart),'k');
    end 
    % meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    title(sprintf('%s normalized ITPC',subjectNames{i})) 
    titleText = sprintf('normalizedITPC_%s',subjectNames{i}); 
    print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
    close all
end

%% plot itpc 
ylims = [1,50]; 
xlims = [size(toi,1),size(toi,2)]; 
toi = p.tstart:1:p.tstop;
ytick = 10:10:numel(foi);
xtick = 1:1000:numel(toi);

for i = 1:numel(subjectNames)
    figure
    set(gcf, 'Position',  [100, 100, 400, 300])
    hold on
    imagesc(A.subject{i})
    xlim(xlims)
    ylim(ylims)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar
    caxis([0,0.6])

    set(gca,'YDir','normal')
     set(gca,'XTick',xtick)
     set(gca,'XTickLabel',toi(xtick))
     set(gca,'YTick',ytick)
     set(gca,'YTickLabel',foi(ytick))
     hold on
    for iEv = 1:numel(p.eventTimes)
        vline(p.eventTimes(iEv)+abs(p.tstart),'k');
    end 
    % meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
    title(sprintf('%s ITPC',subjectNames{i})) 
    titleText = sprintf('ITPC_%s',subjectNames{i}); 
    % print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))
    close all
end

%% plot up down group norm 
ylims = [1,50]; 
xlims = [size(toi,1),size(toi,2)]; 
toi = p.tstart:1:p.tstop;
ytick = 10:10:numel(foi);
xtick = 1:1000:numel(toi);

figure
set(gcf, 'Position',  [100, 100, 400, 300])
hold on
imagesc(A.normDown)
xlim(xlims)
ylim(ylims)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
colorbar
caxis([-0.3,0.3])

set(gca,'YDir','normal')
set(gca,'XTick',xtick)
set(gca,'XTickLabel',toi(xtick))
set(gca,'YTick',ytick)
set(gca,'YTickLabel',foi(ytick))
hold on
for iEv = 1:numel(p.eventTimes)
    vline(p.eventTimes(iEv)+abs(p.tstart),'k');
end
% meg_timeFreqPlotLabels(toi,foi,xtick,ytick,p.eventTimes);
title(sprintf('Downers normalized ITPC')) 
% print(gcf,'-dpng','-painters','normalizedITPC_downers.png')

%% try stats on the image 

im = A.normUp - A.normDown; 
imabs = abs(im); 
thresh = 0; stat = 'mean'; 
[clusterStats, maxAbsClusterStat, C] = rd_clusterStat(imabs, thresh, stat);

%% try stats on the image 

im = A.up - A.down; 

h = zeros(size(im,1),size(im,2));
pval = zeros(size(im,1),size(im,2));
for r = 1:size(im,1)
    for c = 1:size(im,2)
        if all(~isnan(im(r,c))) && all(~isnan(im(r,c)))
            [h(r,c),pval(r,c)] = ttest2(im(r,c),im(r,c));
        else
            h(r,c) = NaN;
            pval(r,c) = NaN;
        end
    end
end

thresh = im>0; 
[clusterStats, maxAbsClusterStat, C] = rd_clusterStat(imabs,thresh);

%% notes from November 11 meeting w Rachel 

parfor i =1:100 % run multicore 
    a = i+1; 
end

%% 

data1 = A.normSubject(:,:,ITPCsubject==1);
data2 = A.normSubject(:,:,ITPCsubject==-1);

%% permutation random shuffling subject level 
% November 11, 2020 

idxUp = find(ITPCsubject==1); 
idxDown = find(ITPCsubject==-1); 
idxAll = [idxUp idxDown]; 

hAll = []; pvalAll = []; 
clusterStats = []; maxAbsClusterStatAll = []; 
for i = 2:1000 % number of random shuffles 
    idxPerm = randperm(size(idxAll,2));
    idxRandUp = idxPerm(1:4);
    idxRandDown = idxPerm(5:8);
    whichGroup = round(rand); % 0 or 1 
    if whichGroup == 1
        idxRandUp = [idxRandUp idxAll(9)];
    else 
        idxRandDown = [idxRandDown idxAll(9)];
    end

%     idxRandUp = idxAll(1:4);
%     idxRandDown = idxAll(5:8);

    upVals = A.normSubject(:,:,idxRandUp); % cat(3,upVals,groupITPCspec.normSubject{idxRandUp(j)});
    downVals = A.normSubject(:,:,idxRandDown); % cat(3,downVals,groupITPCspec.normSubject{idxRandDown(j)});
    
    % point by point ttest (or ranksum for Wilcoxon rank test? equivalent to
    % Mann-Whitney U-test for unequal samples)
    h = zeros(size(upVals,1),size(upVals,2)); 
    pval = zeros(size(upVals,1),size(upVals,2)); 
    for r = 1:size(upVals,1)
        for c = 1:size(upVals,2)
            if all(~isnan(squeeze(upVals(r,c,:)))) && all(~isnan(squeeze(downVals(r,c,:))))
                [h(r,c),pval(r,c)] = ttest2(squeeze(upVals(r,c,:)),squeeze(downVals(r,c,:))); 
            else
                h(r,c) = NaN; 
                pval(r,c) = NaN; 
            end
        end
    end
    
    im = pval;
    thresh = pval<0.05;
    stat = 'sum';
    [clusterStats{i}, maxAbsClusterStat, C(i)] = rd_clusterStat(im, thresh, stat);
    
    % save iteration statistics 
    hAll = cat(3,hAll,h);
    pvalAll = cat(3,pvalAll,pval);
    % clusterStatsAll = cat(1,clusterStatsAll,clusterStats); 
    maxAbsClusterStatAll = cat(1,maxAbsClusterStatAll,maxAbsClusterStat);
    disp(i) % stopped 522 
end

%% plot permuted hypothesis and pvals 

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

titleText = 'ITPCpermutation_clusterStat_histogram_normalized'; % 'ITPCpermutation'; 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% plot cluster distribution 

clusterCat = []; 
for i = 1:500
    clusterCat = [clusterCat clusterStats{i}]; 
end

figure
subplot 211
histogram(clusterCat,100)
cluster95th = prctile(clusterCat,95);
xline(cluster95th,'r',num2str(cluster95th,4))
ylabel('Cluster stat')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off

subplot 212
histogram(maxAbsClusterStatAll,100)
clusterMax95th = prctile(maxAbsClusterStatAll,95);
xline(clusterMax95th,'r',num2str(clusterMax95th,4))
ylabel('Max abs cluster stat')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off

%% use fieldtrip cluster based permutation tests 
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = toi/1000;                  % time window sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg, dataFIC);

TFRhann = x; 

%         label: {149x1 cell}                % Channel names
%        dimord: 'chan_freq_time'            % Dimensions contained in powspctrm, channels X frequencies X time
%          freq: [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30]  % Array of frequencies of interest (the elements of freq may be different from your cfg.foi input depending on your trial length)
%          time: [1x41 double]               % Array of time points considered
%     powspctrm: [149x15x41 double]          % 3-D matrix containing the power values
%          elec: [1x1 struct]                % Electrode positions etc
%          grad: [1x1 struct]                % Gradiometer positions etc
%           cfg: [1x1 struct]                % Settings used in computing this frequency decomposition


%% 
figure
set(gcf, 'Position',  [100, 100, 500, 300])
xtick = [2001,2300+2001]; 
for i = 1 
    hold on
    imagesc(hAll(:,:))
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
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
end

titleText = 'ITPCcluster_TrueData'; % 'ITPCpermutation'; 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))

%% true cluster sizes 
figure
histogram(clusterStats{1},100)
xline(cluster95th,'r',num2str(cluster95th,4))
ylabel('Cluster stat (true data)')
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize = 12;
box off

titleText = 'ITPCcluster_TrueData_histogram'; % 'ITPCpermutation'; 
print(gcf,'-dpng','-painters',sprintf('%s.png',titleText))


