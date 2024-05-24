function meg_manuscriptFigs_normPeak_bar
% function meg_manuscriptFigs_normPeak_bar

%% Load peak data
user = 'kantian';
filename = sprintf('/Users/%s/Dropbox/Data/TANoise/fromRachel/itpcNorm_TS_Peaks_N10_20211225_workspace.mat',user); 
load(filename)

%% Figure settings
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size
showN = 1; % show n = X annotation
figFormat = 'svg'; % svg
annotateMean = 0;
annotateStats = 1;
saveFigs = 1;
showLegend = 1; 
cueLevel = {'cueT1','cueT2'};
plotStyle = 'scatter'; % bar
plotSubjects = 0; 

% Figure directory
[figDir,dateStr,style,colors,p] = meg_manuscriptParams; 

%% Plot bars of normalized ITPC peak by target by precue 
figure
fh = subplot(1,1,1);
hold on
meg_figureStyle
set(gcf,'Position',[100 100 300 style.height])
for iT = 1:2 % target 
    for iC = 1:numel(cueLevel) % precue 
        if iC==1
            x = iT-style.xBufferSml;
        elseif iC==2
            x = iT+style.xBufferSml;
        else
            x = iT; 
        end
        y = peakMean(iC,iT);
        % errorbar(1:2, peakMean(iCue,:), peakDiffSte, '.', 'MarkerSize', 30)

        switch cueLevel{iC}
            case 'all'
                faceColor = colors.lightPurple;
                erColor = colors.darkPurple;
            case 'cueT1'
                faceColor = p.cueColors(iC,:);
                erColor = colors.darkestBlue;
            case 'cueT2'
                faceColor = p.cueColors(iC,:);
                erColor = colors.darkestRed;
        end
        err = peakDiffSte(iT);

        switch plotStyle
            case 'scatter'
                errorbar(x,mean(y,'omitnan'),err,'Marker','.','MarkerSize',style.scatter.MarkerSize,'MarkerFaceColor',faceColor,'MarkerEdgeColor',faceColor,...
                    'Color',faceColor,'LineWidth',2);
                % ylim([0.05 0.1])
            case 'bar'
                pBar = bar(x,mean(y,'omitnan'));
                pBar.BarWidth = 0.34;
                pBars(iC) = pBar;

                er = errorbar(x,y,err,err);
                er.LineWidth = 2;
                er.CapSize = 0;
                er.LineStyle = 'none';

                pBar.EdgeColor = faceColor;
                pBar.FaceColor = faceColor;
                ylim([-0.02 0.1])
        end

        % Plot subjects
        if plotSubjects
            subjectPeak = mean(peakData,4,'omitnan'); 
            yS = squeeze(subjectPeak(iC,iT,:)); 
            yS = yS(~isnan(yS));
            scatter(repmat(x,size(yS)),yS,style.scatter.MarkerSizeS,'filled','MarkerFaceColor','w')
            scatter(repmat(x,size(yS)),yS,style.scatter.MarkerSizeS,'filled','MarkerFaceColor',faceColor,'MarkerEdgeColor','w','MarkerFaceAlpha',0.5)
        end

        % Mean annotation
        if annotateMean
            lbl = sprintf('%0.2f',mean(y,'omitnan'));
            txt = text(x, 0, lbl, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
            txt.Color = [1 1 1];
            txt.FontSize = style.txtSize_Annotation;
            txt.FontName = 'Helvetica-Oblique';
        end

    end
end

% Axes labels
ylabel('Normalized ITPC Peak')

xlabel('Target') 
xticks([1 2])
xticklabels({'T1','T2'})
xlim([1-style.xBuffer/1.5 2+style.xBuffer/1.5])

% Stats annotation
if annotateStats
    txt = meg_annotateStats(1,max(fh.YLim)*0.93,'***'); % T1 
    txt = meg_annotateStats(2.03,max(fh.YLim)*0.95,'ns'); % T2
end

% Legend
if showLegend
    if plotSubjects
        txt = text(max(fh.XLim)*0.98,max(fh.YLim)*1.1,'Precue T1','HorizontalAlignment','right','Color',p.cueColors(1,:),'FontSize',style.txtSize_Legend);
        txt = text(max(fh.XLim)*0.98,max(fh.YLim)*1.01,'Precue T2','HorizontalAlignment','right','Color',p.cueColors(2,:),'FontSize',style.txtSize_Legend);
    else
        txt = text(max(fh.XLim)*0.98,max(fh.YLim)*1.042,'Precue T1','HorizontalAlignment','right','Color',p.cueColors(1,:),'FontSize',style.txtSize_Legend);
        txt = text(max(fh.XLim)*0.98,max(fh.YLim)*1.01,'Precue T2','HorizontalAlignment','right','Color',p.cueColors(2,:),'FontSize',style.txtSize_Legend);
    end
end

if saveFigs
    if plotSubjects
        figTitle = sprintf('meg_manuscriptFigs_normPeak_bar_%s_subjects',dateStr);
    else
        figTitle = sprintf('meg_manuscriptFigs_normPeak_bar_%s',dateStr);
    end
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end



