% paper_figures_2019.m

%% setup
gd = '/Users/rachel/Google Drive';
% gd = '/Local/Users/denison/Google Drive';

dataDir = sprintf('%s/NYU/Projects/Temporal_Attention/Code/Expt_Scripts/Behav/data', gd);
fitDir = sprintf('%s/NYU/Projects/Temporal_Attention/Code/Models/TA_Model/fit/startpoints', gd);
mrDir = sprintf('%s/NYU/Projects/Temporal_Attention/Code/Models/TA_Model/fit/modelrecovery', gd);
figDir = sprintf('%s/NYU/Manuscripts/TA_Model/Figures', gd);

dataFile = 'E2_SOA_cbD6_run98_N5_norm_workspace_20191003.mat';
D = load(sprintf('%s/%s', dataDir, dataFile));

% DD = load([dataDir '/E2_SOA_cbD6_run98_N5_workspace_20180731.mat'])

%% more setup
modelClass0 = 'valo';
switch modelClass0
    case 'valo'
        spDir1 = '20191010_valo_12params';
        spDir2 = '20191016_valo_nolimit';
        
        fitPath1 = sprintf('%s/%s/*-34_valo.mat', fitDir, spDir1);
        fitPath2 = sprintf('%s/%s/*-4_valo_nolimit.mat', fitDir, spDir2);
        
        resampFile1 = 'fit/resample/20191005_valo/resampAnalysis_20191103.mat';
        resampFile2 = 'fit/resample/20191105_valo_nolimit/resampAnalysis_20191112.mat';
    case 'span'
        spDir1 = '20191010_span_15params';
        spDir2 = '20191016_span_nolimit';
        
        fitPath1 = sprintf('%s/%s/*-37_span.mat', fitDir, spDir1);
        fitPath2 = sprintf('%s/%s/*-40_span_nolimit.mat', fitDir, spDir2);
    case '1-attK'
        spDir1 = '20191015_1-attK_14params';
        spDir2 = '20191016_1-attK_nolimit';
        
        fitPath1 = sprintf('%s/%s/*-6_1-attK.mat', fitDir, spDir1);
        fitPath2 = sprintf('%s/%s/*-22_1-attK_nolimit.mat', fitDir, spDir2);
    case '3S'
        spDir1 = '20191015_3S_12params';
        spDir2 = '20191016_3S_nolimit';
        
        fitPath1 = sprintf('%s/%s/*-26_3S.mat', fitDir, spDir1);
        fitPath2 = sprintf('%s/%s/*-26_3S_nolimit.mat', fitDir, spDir2);
    case 'valov'
        spDir1 = '20191102_valov_9params';
        spDir2 = '20191105_valov_nolimit';
        
        fitPath1 = sprintf('%s/%s/*-35_valo.mat', fitDir, spDir1);
        fitPath2 = sprintf('%s/%s/*-36_valov_nolimit.mat', fitDir, spDir2);
    otherwise
        error('modelClassName not recognized')
end

plotTaskArray = false;

if strcmp(modelClass0,'3S')
    x = 1;
else
    x = 0;
end

%% load fits
f1 = dir(fitPath1);
f2 = dir(fitPath2);

fitFile1 = sprintf('%s/%s/%s',fitDir, spDir1, f1.name);
fitFile2 = sprintf('%s/%s/%s',fitDir, spDir2, f2.name);

fit1 = load(fitFile1);
fit2 = load(fitFile2);

%% task array
if plotTaskArray
rsoa = [3 9 10];
nsoa = numel(rsoa);

rcond = 2:4;
ncond = numel(rcond);

% set params
opt.neutralT1Weight = .5;
opt.span = 800;

for isoa = 1:nsoa
    for icond = 1:ncond
        [perfv, p(isoa,icond), ev] = runModelTA([], [], rsoa(isoa), [], rcond(icond));
    end
end

h = figure('Position',[500 80 750 230]);
for isoa = 1:nsoa
    for icond = 1:ncond
        task = p(isoa,icond).task(1,:);
        
        subplot(ncond,nsoa,(icond-1)*nsoa+isoa)
        plot(p(1,1).tlist, task)
    end
end

ax = h.Children;

xlims = [300 1700];
ylims = [0 1];

% set axis properties
for i = 1:numel(ax)
    ax(i).XLim = xlims;
    ax(i).YLim = ylims;
    ax(i).Visible = 'off';
%     ax(i).FontSize = 14;
%     ax(i).LineWidth = 1;
end
for i = 1:3
    ax(i).Children(1).Color = [0 0 0];
end
for i = 4:6
    ax(i).Children(1).Color = [.27 .27 .27];
end
for i = 7:9
    ax(i).Children(1).Color = [.4 .4 .4];
end
ax(7).Visible = 'on';
ax(7).Box = 'off';
ax(7).XTick = [];
ax(7).YTick = [0 1];
ax(7).LineWidth = 1;
ax(7).TickDir = 'out';

% print_pdf('task_t1t2n_soa200-500-800.pdf', figDir)
end

%% time series: model 1
% turn plotTS on in runModel

% set params
modelClass = fit1.modelClass;
opt = fit1.opt;
% opt.neutralT1Weight = 0.5;
% opt.vAttScale2 = 1;

rsoa = 5; % 300
rcond = 2; % endoT1

[perfv, p, ev] = runModel(opt, modelClass, rsoa, [], rcond);

% avAmp = max(p.rav(:))*p.aAV;
% aiAmp = max(p.rai(:))*p.aAI;
% avaiRatio = avAmp/aiAmp;

h = gcf;
ax = h.Children;
s = 1.3; % y-axis scale factor

ax(6+x).YLim(1) = min(p.rav(:))*s; % AV
ax(5+x).YLim(1) = min(p.rai(:))*s; % AI
ax(4+x).YLim(1) = min(p.r1(:))*s; % S1
ax(3+x).YLim(1) = min(p.r2(:))*s; % S2
ax(2+x).YLim(1) = min(-p.rd(:))*2; % D

ax(7+x).YLim(2) = max(p.stim(:))*s; % stim
ax(6+x).YLim(2) = max(p.rav(:))*s; % AV
ax(5+x).YLim(2) = max(p.rai(:))*s; % AI
ax(4+x).YLim(2) = max(p.r1(:))*s; % S1
ax(3+x).YLim(2) = max(p.r2(:))*s; % S2
ax(2+x).YLim(2) = max(-p.rd(:))*2; % D

ax(6+x).YLabel.String = {'Voluntary','attention'}; 
ax(5+x).YLabel.String = {'Involuntary','attention'};
ax(4+x).YLabel.String = 'Sensory 1';
ax(3+x).YLabel.String = 'Sensory 2';

if strcmp(modelClass,'3S')
    ax(2+x).YLim(1) = min(p.r3(:))*s; % S3
    ax(2+x).YLim(2) = max(p.r3(:))*s; % S3
    ax(2+x).YLabel.String = 'Sensory 3';
end

ax(6+x).Children(1).Color = [71 98 159]/255;
ax(5+x).Children(1).Color = [112 153 208]/255;

for i = 1:numel(ax(3+x).Children)
    ax(3+x).Children(i).Color = [34 139 34]/255; % med green
end
if strcmp(modelClass,'3S')
    for i = 1:numel(ax(3).Children)
        ax(3).Children(i).Color = [59 117 56]/255; % dark green
    end
end
for i = 1:numel(ax(2).Children)
    ax(2).Children(i).Color = [110 110 110]/255;
end
% flip decision variable to positive
ax(2).Children(3).YData = -ax(2).Children(3).YData;
ax(2).Children(4).YData = -ax(2).Children(4).YData;

h.Children(1).Title.String = ''; % turn off tile
h.Position(3) = 330; % make wider
for i = 1:numel(ax)
    ax(i).TickDir = 'out';
    ax(i).LineWidth = 0.5;
    ax(i).YTick = [];
    ax(i).YLabel.Rotation = 0;
    ax(i).YLabel.HorizontalAlignment = 'right';
    ax(i).YLabel.VerticalAlignment = 'middle';
end

% print_pdf(sprintf('ts_%s%s_soa300_endoT1_%s', ...
%     fit1.modelClass, fit1.jobStr, datestr(now,'yyyymmdd')), figDir)

for i = 2:numel(ax)
    ax(i).FontSize = 8.3;
    ax(i).LabelFontSizeMultiplier = 9/8.3;
end

aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [6 6*aspectRatio];
h.PaperPosition = [0 0 6 6*aspectRatio];
print(sprintf('%s/ts_%s%s_soa300_endoT1_%s', figDir, ...
    fit1.modelClass, fit1.jobStr, datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% dprime data only
h = plotFit(D, fitFile1, [], 1, 1, [], 'paper');

% print_pdf(sprintf('plot_dprime_%s', datestr(now,'yyyymmdd')), figDir)

ax = h.Children;
for i = [1 2 4 5]
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
ax(3).FontSize = 8.3; % legend

ax(1).YTickLabel{1} = '0.5';
ax(2).YTickLabel{1} = '0.5';

aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [8.5 8.5*aspectRatio];
h.PaperPosition = [0 0 8.3 8.3*aspectRatio];
print(sprintf('%s/plot_dprime_%s', figDir, datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% dprime cueing effect scatterplots for each SOA
% use unnormalized data?
% dataFile = 'E2_SOA_cbD6_run98_N5_workspace_20191003.mat';
targetColors = {[137 65 143]/255,[0 165 133]/255};

measure = 'RT';
switch measure
    case 'dprime'
        vals = D.dpData;
        lims = [-.5 3.5];
        xlab = 'Invalid performance (d'')';
        ylab = 'Valid performance (d'')';
    case 'RT'
        vals = D.rtData;
        lims = [.5 1.11];
        xlab = 'Invalid RT (s)';
        ylab = 'Valid RT (s)';
end

h = figure('Position',[200 200 1000 400]);
for iSOA = 1:10
    subplot(2,5,iSOA)
    hold on
    plot(lims,lims,'k','LineWidth',.75)
    for iT = 1:2
        valsV = squeeze(vals{iT}(1,iSOA,:));
        valsI = squeeze(vals{iT}(2,iSOA,:));
        p1(iT) = plot(valsI, valsV, '.','MarkerSize',8,'Color',targetColors{iT});
    end
    xlim(lims)
    ylim(lims)
    axis square
    title(sprintf('%d', D.t1t2soa(iSOA)))
    set(gca,'XTick',[.5 .7 .9 1.1])
    set(gca,'YTick',[.5 .7 .9 1.1])
    if iSOA==6
        xlabel(xlab)
        ylabel(ylab)
    else
        set(gca,'XTickLabels',[])
        set(gca,'YTickLabels',[])
    end
end
legend(p1,{'T1','T2'})
legend boxoff

ax = h.Children;
ax(1).FontSize = 8.3; % legend
ax(1).Position(1) = .9;
for i = 2:numel(ax)
    ax(i).LineWidth = 0.5;
    ax(i).FontSize = 8.3;
    ax(i).TickDir = 'out';
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [17 17*aspectRatio];
h.PaperPosition = [0 0 17 17*aspectRatio];
print(sprintf('%s/plot_%sSOAScatterVINorm_%s', figDir, measure, datestr(now,'yyyymmdd')),'-dpdf','-r0')
    
%% dprime benefit and cost scatterplots for each SOA and target
% use unnormalized data?
% dataFile = 'E2_SOA_cbD6_run98_N5_workspace_20191003.mat';
lims = [-.5 3.5];
for iT = 1:2
    h = figure('Position',[200 200 1000 400]);
    for iSOA = 1:10
        subplot(2,5,iSOA)
        hold on
        plot(lims,lims,'k','LineWidth',.75)
        valsV = squeeze(D.dpData{iT}(1,iSOA,:));
        valsI = squeeze(D.dpData{iT}(2,iSOA,:));
        valsN = squeeze(D.dpData{iT}(3,iSOA,:));
        p1(1) = plot(valsN, valsV, '.','MarkerSize',20);
        p1(2) = plot(valsN, valsI, '.','MarkerSize',20);
        xlim(lims)
        ylim(lims)
        axis square
        title(sprintf('%d', D.t1t2soa(iSOA)))
        if iSOA==6
            xlabel('d'' neutral')
            ylabel('d'' valid or invalid')
        else
            set(gca,'XTickLabels',[])
            set(gca,'YTickLabels',[])
        end
    end
    legend(p1,{'V vs. N','I vs. N'})
    legend boxoff
end

ax = h.Children;
ax(1).FontSize = 8.3; % legend
ax(1).Position(1) = .9;
for i = 2:numel(ax)
    ax(i).LineWidth = 0.5;
    ax(i).FontSize = 8.3;
    ax(i).TickDir = 'out';
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [17 17*aspectRatio];
h.PaperPosition = [0 0 17 17*aspectRatio];
% print(sprintf('%s/plot_dprimeSOAScatterVINorm_%s', figDir, datestr(now,'yyyymmdd')),'-dpdf','-r0')
    
%% RT
rtFields = {'rtMean','rtSte','rtDataCueEff','rtData'};
f1 = plotFit(D, fitFile1, [], 1, 1, rtFields, 'paper');

f2 = figure('Position',[100 100 580 250]);
h = copyobj(f1.Children(3:5), f2);
h(3).YLabel.String = 'RT (s)';
h(3).XLabel.String = 'SOA (ms)';
for i = 2:3
    h(i).YLim = [.65 .95];
end

% print_pdf(sprintf('plot_RT_%s', datestr(now,'yyyymmdd')), figDir)

ax = h;
for i = [2 3]
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
ax(1).Visible = 'off'; % legend

ax(2).YTick = .7:.1:.9;
ax(3).YTick = .7:.1:.9;

h = gcf;
aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [8.5 8.5*aspectRatio];
h.PaperPosition = [0 0 8.3 8.3*aspectRatio];
print(sprintf('%s/plot_RT_%s', figDir, datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% fit: model1
plotFit(D, fitFile1, [], [], [], [], 'paper')
% plotFit(D, fitFile1, resampFile1, [], [], [], 'paper')

f1 = gcf;
f1.Children(5).Children(end-1).LineStyle = '--'; % T1 invalid
f1.Children(4).Children(end-2).LineStyle = '--'; % T2 valid

% print_pdf(sprintf('plot_%s%s_%s', ...
%     fit1.modelClass, fit1.jobStr, datestr(now,'yyyymmdd')), figDir)

h = f1;
ax = h.Children;
for i = [1 2 4 5]
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
ax(3).FontSize = 8.3; % legend

ax(1).YTickLabel{1} = '0.5';
ax(2).YTickLabel{1} = '0.5';

aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [8.5 8.5*aspectRatio];
h.PaperPosition = [0 0 8.3 8.3*aspectRatio];
print(sprintf('%s/plot_%s%s_%s', figDir, fit1.modelClass, fit1.jobStr, ...
    datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% fit: model2
% plotFit(D, fitFile2, [], [], [], [], 'paper')
plotFit(D, fitFile2, resampFile2, [], [], [], 'paper')

f2 = gcf;
f2.Children(5).Children(end-2).LineStyle = '--'; % T1 valid
f2.Children(4).Children(end-2).LineStyle = '--'; % T2 valid

% print_pdf(sprintf('plot_%s%s_%s', ...
%     fit2.modelClass, fit2.jobStr, datestr(now,'yyyymmdd')), figDir)

h = f2;
ax = h.Children;
for i = [1 2 4 5]
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
ax(3).FontSize = 8.3; % legend

ax(1).YTickLabel{1} = '0.5';
ax(2).YTickLabel{1} = '0.5';

aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [8.5 8.5*aspectRatio];
h.PaperPosition = [0 0 8.3 8.3*aspectRatio];
print(sprintf('%s/plot_%s%s_%s', figDir, fit2.modelClass, fit2.jobStr, ...
    datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% calculate AIC
% total params, which can be more than the number of params optimized
switch modelClass0
    case 'valo'
        nParams = 20;
    case 'span'
        nParams = 23;
    case '1-attK'
        nParams = 21;
    case '3S'
        nParams = 22;
    case 'valov'
        nParams = 16;
    otherwise
        error('modelClass not found')
end

nParamsNL = nParams - 2;

n = 60; % number of data points
% fit = load(fitFile1, 'opt','R2','cost');
% k = numel(fields(fit1.opt)); % number of parameters
k = nParams;
aic = 2*k + n*log(fit1.cost);
fprintf('%s, %s, %d params: R2 = %.04f, AIC - constant = %.2f\n\n', ...
    fit1.modelClass, fit1.jobStr(2:end), nParams, fit1.R2, aic)

% fit = load(fitFile2, 'opt','R2','cost');
% k = numel(fields(fit2.opt)); % number of parameters
k = nParamsNL; % assume fit2 is the nolimit version of fit1
aic = 2*k + n*log(fit2.cost);
fprintf('%s, %s, %d params: R2 = %.04f, AIC - constant = %.2f\n\n', ... 
    fit2.modelClass, fit2.jobStr(2:end), nParamsNL, fit2.R2, aic)


%% fit of Denison 2017 Exp 1
% colors
colors = get(gcf,'DefaultAxesColorOrder');
colors(3,:) = [.5 .5 .5];
c = rgb2hsv(colors);
c(1:2,2) = .2;
c(1:2,3) = .95;
c(3,3) = .9;
barcolors = hsv2rgb(c);
c(1:2,2) = .4;
c(1:2,3) = .95;
c(3,3) = .7;
ebcolors = hsv2rgb(c);

fitDir0 = sprintf('%s/NYU/Projects/Temporal_Attention/Code/Models/TA_Model/fit', gd);

% fitName0 = 'fit_workspace_20200807T1633_valo.mat'; % 3 params
fitName0 = 'fit_workspace_20200807T1719_valo.mat'; % 2 params
% fitName0 = 'fit_workspace_20200808T1017_valo.mat'; % 3 params, averageContrasts
fitFile0 = sprintf('%s/%s',fitDir0, fitName0);

fit0 = load(fitFile0);

dataFile0 = 'denison2017_exp1';
% dataFile0 = 'denison2017_exp1_averageContrasts';
D0 = load(sprintf('%s/%s.mat', dataDir, dataFile0));
D0 = D0.D;

% plot
targetNames = {'T1','T2'};
cueNames = {'V','I','N'};
dpLims = [0 2.5];
idx = [1 3 2];
nT = 2;
nV = 3;

figure
for iT = 1:nT
    subplot(1,2,iT);
    hold on
    
    for iV = 1:nV
        i = idx(iV);
        b1(iV) = bar(iV, D0.dpMean{iT}(i),'FaceColor',barcolors(i,:),...
            'EdgeColor','none');
        p1(iV) = errbar(iV, D0.dpMean{iT}(i), D0.dpSte{iT}(i),...
            'LineWidth',1.5,'Color',ebcolors(i,:));
        p2(iV) = plot(iV, fit0.model(i,1,iT), 'o', 'MarkerSize',8,...
            'Color',colors(i,:));
    end
    
    set(gca,'XTick',1:nV)
    set(gca,'XTickLabel', cueNames(idx))
    ylim(dpLims)
    set(gca,'TickDir','out','LineWidth',0.5)
    
    if iT==1
        xlabel('Precue')
        ylabel('d''')
    end
    title(targetNames{iT})
    box off
end

print(sprintf('%s/plot_%s_%s%s_%s', figDir, dataFile0, fit0.modelClass, fit0.jobStr, ...
    datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% fit of Fernandez 2019
fitDir0 = sprintf('%s/NYU/Projects/Temporal_Attention/Code/Models/TA_Model/fit', gd);

visFieldLoc = 'all';
switch visFieldLoc
    case 'fov'
        fitName0 = 'fit_workspace_20200807T1738_valo.mat'; % 2 params   
    case 'rh'
        fitName0 = 'fit_workspace_20200807T1749_valo.mat'; % 2 params
    case 'uvm'
        fitName0 = 'fit_workspace_20200807T1758_valo.mat'; % 2 params
    case 'all'
        fitName0 = 'fit_workspace_20200808T1030_valo.mat'; % 2 params
end
fitFile0 = sprintf('%s/%s',fitDir0, fitName0);

fit0 = load(fitFile0);

dataFile0 = sprintf('fernandez2019_%s', visFieldLoc');
D0 = load(sprintf('%s/%s.mat', dataDir, dataFile0));
D0 = D0.D;

% plot
targetNames = {'T1','T2'};
cueNames = {'V','I','N'};
dpLims = [0 2.5];
idx = [1 3 2]; % [1 3]
nT = 2;
nV = 3; % 2

figure
for iT = 1:nT
    subplot(1,2,iT);
    hold on
    
    for iV = 1:nV
        i = idx(iV);
        b1(iV) = bar(iV, D0.dpMean{iT}(i),'FaceColor',barcolors(i,:),...
            'EdgeColor','none');
        p1(iV) = errbar(iV, D0.dpMean{iT}(i), D0.dpSte{iT}(i),...
            'LineWidth',1.5,'Color',ebcolors(i,:));
        p2(iV) = plot(iV, fit0.model(i,1,iT), 'o', 'MarkerSize',8,...
            'Color',colors(i,:));
    end
    
%     b1 = bar(1:nV, D0.dpMean{iT}(idx),'FaceColor',[.7 .7 .7]);
%     p1 = errbar(1:nV, D0.dpMean{iT}(idx), D0.dpSte{iT}(idx),...
%         'LineWidth',1.5,'Color',[.3 .3 .3]);
%     p2 = plot(1:nV, fit0.model(idx,1,iT), 'o', 'Color','k');
    
    set(gca,'XTick',1:nV)
    set(gca,'XTickLabel', cueNames(idx))
    ylim(dpLims)
    set(gca,'TickDir','out','LineWidth',0.5)
    
    if iT==1
        xlabel('Precue')
        ylabel('d''')
    end
    title(targetNames{iT})
    box off
end

print(sprintf('%s/plot_%s_%s%s_%s', figDir, dataFile0, fit0.modelClass, fit0.jobStr, ...
    datestr(now,'yyyymmdd')),'-dpdf','-r0')

%% model recovery
mrFile1 = sprintf('%s/20200808_MR_gen-valo_workspace.mat', mrDir);
mrFile2 = sprintf('%s/20200808_MR_gen-valo_nolimit_workspace.mat', mrDir);

mr1 = load(mrFile1);
mr2 = load(mrFile2);

relAIC1 = mr1.sortedAIC - min(mr1.sortedAIC(:));
relAIC2 = mr2.sortedAIC - min(mr2.sortedAIC(:));

fitNames = {'Main model fit', 'No limit model fit'};

colors = [0 0 0; .7 .7 .7];

h = figure;
ax(1) = subplot(1,2,1);
ax(1).ColorOrder = colors; hold on
plot(relAIC1, 'LineWidth',1)
xlabel('Fit number')
ylabel('\DeltaAIC')
title({'Main model','generated'})
axis square
ax(2) = subplot(1,2,2);
ax(2).ColorOrder = colors; hold on
plot(relAIC2, 'LineWidth',1)
title({'No limit model','generated'})
legend(fitNames,'Location','SE')
legend boxoff
axis square

for i = 1:numel(ax)
    ax(i).TickDir = 'out';
    ax(i).FontSize = 8.3;
    ax(i).LabelFontSizeMultiplier = 9/8.3;
    ax(i).TitleFontSizeMultiplier = 9/8.3;
end
h.Children(1).FontSize = 8.3; % legend

aspectRatio = h.Position(4)/h.Position(3);
h.PaperUnits = 'centimeters';
h.PaperSize = [8.5 8.5*aspectRatio];
h.PaperPosition = [0 0 8.3 8.3*aspectRatio];
print(sprintf('%s/plot_modelRecovery_%s', figDir, ...
    datestr(now,'yyyymmdd')),'-dpdf','-r0')

