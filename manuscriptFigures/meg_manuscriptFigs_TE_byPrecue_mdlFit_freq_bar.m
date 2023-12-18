function meg_manuscriptFigs_TE_byPrecue_mdlFit_freq_bar
% function meg_manuscriptFigs_TE_byPrecue_mdlFit_freq_bar
% Plot model fits (linear+periodic) free frequency by precue (All trials,
% precue T1, precue T2) 

%% Load data 
user = 'kantian'; 
filename = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ModelFit_1_freeFreq_231103.mat',user);
load(filename)

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotSubjects = 0; 
annotateMean = 0; 
annotateStats = 1; 
saveFigs = 1; 
cueLevel = {'all','cueT1','cueT2'}; 
[style, colors] = meg_manuscriptStyle;

% Figure directory 
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/manuscriptFigures/figs',user); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Find solution w minimum fval 
for iF = 2 % just linear+periodic model 
    for iS = 1:size(data.all.(fitLevel),2) % sessions
        for iC = 1:numel(cueLevel)
            % --- Plot model fit ---
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS)  );
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:) = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:);
        end
    end
end

%% Plot bars of fitted frequency
% --- Plot figure (bar): Fitted frequency by precue  --- 
figure
fh = subplot(1,1,1); 
hold on
meg_figureStyle
set(gcf,'Position',[100 100 200 style.height])
for iC = 1:numel(cueLevel)
    clear x y idx paramNames
    paramNames = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames;
    paramOI = 'freq';
    idx = find(contains(paramNames,paramOI));
    y = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(:,idx); % sessions
    y = meg_sessions2subjects(y'); % subjects 
    x = iC; 

    pBar = bar(x,mean(y,'omitnan'));
    pBar.BarWidth = 0.85;

    % standard error of mean
    err = std(y)/sqrt(numel(y)); 

    er = errorbar(x,mean(y,'omitnan'),err,err);
        er.LineWidth = 2;
        er.CapSize = 0;
        er.LineStyle = 'none';

    darkMode = 1; % ight or dark color palette
    switch cueLevel{iC}
        case 'all'
            if darkMode
                pBar.FaceColor = colors.darkPurple;
                er.Color = colors.lightPurple;
                sColor = colors.mediumPurple; 
            else
                pBar.FaceColor = colors.lightPurple;
                er.Color = colors.darkPurple;
            end
        case 'cueT1'
            if darkMode
                pBar.FaceColor = colors.darkBlue;
                er.Color = colors.lightBlue;
                sColor = colors.precueBlue; 
            else
                pBar.FaceColor = colors.lightBlue;
                er.Color = colors.darkestBlue;
            end
        case 'cueT2'
            if darkMode
                pBar.FaceColor = colors.darkRed;
                er.Color = colors.lightRed;
                sColor = colors.precueRed; 
            else
                pBar.FaceColor = colors.lightRed;
                er.Color = colors.darkestRed;
            end
    end
    pBar.EdgeColor = pBar.FaceColor;

    % Plot subjects scatter points
    if plotSubjects
        scatter(x,y,'filled','MarkerFaceColor',sColor)
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

% Axes labels
ylabel('Frequency (Hz)')

xlabel('Precue') % need to retain sizing?
xticks([1 2 3])
xticklabels({'All','T1','T2'})
xlim([1-style.xBuffer/1.5 3+style.xBuffer/1.5])

% Stats annotation
if annotateStats
    txt = text(2,max(fh.YLim),'n.s.','EdgeColor','none',...
        'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','Bottom','FontName','Helvetica');
    % txt = text(2,max(fh.YLim)*0.95,'**','EdgeColor','none',...
    %     'FontSize',30,'HorizontalAlignment','center','VerticalAlignment','Bottom','FontName','Times');
end

if saveFigs
    figTitle = sprintf('meg_manuscriptFigs_TE_byPrecue_mdlFit_freq_bar_%s',dateStr);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end






