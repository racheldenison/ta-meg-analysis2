function meg_manuscriptFigs_ITPCbyITI_slopes_bar
% function meg_manuscriptFigs_ITPCbyITI_slopes_bar
% load groupITPC_ITI
% and run kt_ITPCbyITI for the model fits by jitter 

%% Settings 
plotSubjects = 0; 
saveFigs = 1; 
annotateStats = 1; 
figFormat = 'svg'; 

% Figure directory
[figDir,dateStr,style,colors,p] = meg_manuscriptParams; 

%% --- Figure (fitted slopes bars by jitter) ---
figure
set(gcf,'Position',[100 100 300 style.height])

fh = subplot(1,1,1);
hold on 
meg_figureStyle

ITIs = 500:200:1500; 
for i = 1:numel(ITIs)
    fieldname = sprintf('ITI%d',ITIs(i));
    y = A.(fieldname).slopes' * 1000;
end

ylim([0 0.1])

xITI = 1:6;
for i = 1:numel(ITIs)
    fieldname = sprintf('ITI%d',ITIs(i));
    y = A.(fieldname).slopes' * 1000;
    y = meg_sessions2subjects(y'); 

    ste = std(y)/sqrt(10);
    meany = mean(y);

    if plotSubjects
        x = repmat(xITI(i),[20,1]);
        scatter(x,y,'filled','MarkerFaceColor',p.cueColors(i,:),'MarkerFaceAlpha',0.5)
    end
    errorbar(xITI(i),meany,ste,'Marker','.','MarkerSize',style.scatter.MarkerSize,'MarkerFaceColor',p.cueColors(i,:),'MarkerEdgeColor',p.cueColors(i,:),...
        'Color',p.cueColors(i,:),'LineWidth',2);


        % pBar = bar(xITI(i),mean(y,'omitnan'));
        % pBar.BarWidth = 0.8;
        % 
        % % standard error of mean
        % err = std(y)/sqrt(numel(y));
        % 
        % er = errorbar(xITI(i),mean(y,'omitnan'),err,err);
        % er.LineWidth = 2;
        % er.CapSize = 0;
        % er.LineStyle = 'none'; 
end
xlim([0 7])
xticks(xITI)
xticklabels({'500','700','900','1100','1300','1500'})
ylabel(sprintf('Slope (\\Delta ITPC/s)'))
xlabel('Jitter (ms)')

% Stats annotation
if annotateStats
    txt = meg_annotateStats( (max(fh.XLim)-min(fh.XLim))/2, max(fh.YLim),'ns');
end

if saveFigs
    figTitle = sprintf('meg_manuscriptFigs_ITPCbyITI_slopes_bar_%s',dateStr);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end
