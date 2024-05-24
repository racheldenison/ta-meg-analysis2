function meg_manuscriptFigs_TE_byPrecue_mdlFit_2Hz_bars
% function meg_manuscriptFigs_TE_byPrecue_mdlFit_2Hz_bars
% Plot 2 Hz model fit 

%% Load data
user = 'kantian'; 
% Model fit data 
filename = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ITPCfit_separate1_2Hz/231213/ModelFit_SeparatePrecueT1T2_2Hz_231213.mat',user);
% ITPC data 
load(filename)

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotErrorBars = 1; % turn off for subject-level plots 
restrictYLim = 1; % turn on for group-level manuscript matching ylims 
plotPrePrecue = 0; % bar showing pre-precue baseline period
plotSubjects = 0; 
saveFigs = 1; 
annotateMean = 0; 
annotateStats = 1; 
plotStyle = 'scatter'; % bar or scatter 

% Figure directory
[figDir,dateStr,style,colors,p] = meg_manuscriptParams; 

% Data settings
cueLevel = {'cueT1','cueT2'}; 
fitLevel = 'session'; 
paramLabelNames = {'Intercept (ITPC)','Slope ','Amplitude (ITPC)','Phase'}; 
paramLabelNames{2} = sprintf('Slope (\\Delta ITPC/s)'); 

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

%% Plot bars of fitted parameter
% --- Plot figure (bar): Fitted parameter by precue  ---
for iP = 1:3
    figure
    fh = subplot(1,1,1);
    hold on
    meg_figureStyle
    set(gcf,'Position',[100 100 150 style.height])
    fh.InnerPosition = [0.4140 0.1570 0.4910 0.7680]; 
    for iC = 1:numel(cueLevel)
        clear x y idx paramNames
        paramNames = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames;
        paramOI = paramNames{iP};
        % Axes labels
        ylabel(paramLabelNames{iP})

        idx = find(contains(paramNames,paramOI));
        y = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(:,idx); % sessions
        y = meg_sessions2subjects(y'); % subjects
        x = iC;

        % standard error of mean
        % err = std(y)/sqrt(numel(y));

        % standard error of difference 
        y1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).minSolution(:,idx); 
        y1 = meg_sessions2subjects(y1'); % subjects
        y2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).minSolution(:,idx); 
        y2 = meg_sessions2subjects(y2'); % subjects
        err = std(y1-y2)/sqrt(numel(y)); 

        [faceColor,erColor,sColor] = meg_manuscriptStyleCue(cueLevel{iC});

        switch plotStyle
            case 'bar'
                pBar = bar(x,mean(y,'omitnan'));
                pBar.BarWidth = 0.85;
                er = errorbar(x,mean(y,'omitnan'),err,err);
                er.LineWidth = 2;
                er.CapSize = 0;
                er.LineStyle = 'none';
                er.Color = erColor;
                pBar.FaceColor = faceColor;
                pBar.EdgeColor = pBar.FaceColor;
            case 'scatter'
                e = errorbar(x,mean(y,'omitnan'),err,'Marker','.','MarkerSize',style.scatter.MarkerSize,'MarkerFaceColor',faceColor,'MarkerEdgeColor',faceColor,...
                    'Color',faceColor,'LineWidth',2);
                e.CapSize = style.scatter.errCapSize; 
        end

        % Plot subjects scatter points
        if plotSubjects
            scatter(x,y,style.scatter.MarkerSizeS,'filled','MarkerFaceColor','w')
            scatter(x,y,style.scatter.MarkerSizeS,'filled','MarkerFaceColor',faceColor,'MarkerEdgeColor','w','MarkerFaceAlpha',0.5)
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
    
    if iP==1
        xlabel('Precue') % need to retain sizing?
    else 
        xlabel('')
    end
    xticks([1 2])
    xticklabels({'T1','T2'})
    xlim([1-style.xBuffer/1.5 2+style.xBuffer/1.5])

    % y axis styling
    if plotSubjects
        switch paramOI
            case 'amplitude'
                yticks(0:0.01:0.05)
                ylim([0 0.05])
        end
    else
        switch paramOI
            case 'amplitude'
                yticks(0.0:0.01:0.04)
                ylim([0.0 0.035])
            case 'slope'
                yticks(0:0.02:0.085)
                ylim([0 0.085])
            case 'intercept'
                yticks(0.1:0.02:0.45)
                ylim([0.26 0.34])
        end
    end

    % Subject lines
    if plotSubjects
        for i = 1:numel(subjectNames)
            s = plot([1 2],[y1(i) y2(i)],'Color',colors.lightgrey);
        end
    end

    % Stats annotation
    if annotateStats
        txt = meg_annotateStats(1.53,max(fh.YLim),'ns');
    end

    if saveFigs
        if plotSubjects
            figTitle = sprintf('meg_manuscriptFigs_TE_byPrecue_mdlFit_2Hz_bar_%s_%s_subjects',paramOI,dateStr);
        else
            figTitle = sprintf('meg_manuscriptFigs_TE_byPrecue_mdlFit_2Hz_bar_%s_%s',paramOI,dateStr);
        end
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
    end
end

%% Export csv long format data for R 
exportCSV = 1;
cueLevel = {'cueT1','cueT2'}; 
clear V
tableHeaders = {'subject','session','precue',...
    'intercept','slope','amplitude'};
varNames = {'intercept','slope','amplitude'}; 
count = 1;
for iS = 1:10 % subjects
    for iSession = 1:2
        for iC = 1:numel(cueLevel)
            for iV = 1:numel(varNames)
                V.subject(count) = iS;
                V.session(count) = iSession;
                V.precue(count)  = iC;
                if iSession==1
                    idxS = iS*2-1; 
                elseif iSession==2
                    idxS = iS*2; 
                end
                V.(tableHeaders{iV+3})(count) = mdlFit.linear2Hz.(cueLevel{iC}).session.minSolution(idxS,iV);
            end
            count = count+1;
        end
    end
end

% --- Export to csv for R ANOVA ---
if exportCSV
    sz = [numel(V.(tableHeaders{1})) numel(tableHeaders)];

    varTypes = repmat("double",[1 numel(tableHeaders)]);
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

    for iT = 1:numel(tableHeaders)
        T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
    end
    
    csvDir = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data'; 
    csvName = sprintf('TANoise_MdlFit_Linear2Hz');
    csvPath = sprintf('%s/%s.csv', csvDir, csvName);
    writetable(T,csvPath)
end






