function meg_manuscriptFigs_TE_allTrials_mdlFit_slope_bar
% meg_manuscriptFigs_TE_allTrials_mdlFit_slope_bar
% Plot linear model slope (All trials) 

%% Load data 
user = 'kantian'; 
filename = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ModelFit_SeparatePrecueT1T2_2Hz_session_231220.mat',user); 
load(filename)

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotSubjects = 0; 
annotateMean = 0; 
annotateStats = 1; 
saveFigs = 1; 
cueLevel = {'all'}; 
plotStyle = 'scatter'; % bar or scatter 

% Figure directory
[figDir,dateStr,style,colors,p] = meg_manuscriptParams; 

%% Find solution w minimum fval 
for iF = 1 % only linear model 
    for iS = 1:size(mdlFit.linear.all.session.solution,2) % sessions
        for iC = 1:numel(cueLevel)
            % --- Plot model fit ---
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS)  );
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:) = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:);
        end
    end
end

%% Plot bar of fitted slope (linear only model) 
figure
fh = subplot(1,1,1); 
hold on
meg_figureStyle
set(gcf,'Position',[100 100 130 style.height])
iF = 1; % linear model 
for iC = 1:numel(cueLevel)
    clear x y idx paramNames
    paramNames = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames;
    paramOI = 'slope';
    idx = find(contains(paramNames,paramOI));
    y = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(:,idx); % sessions
    y = meg_sessions2subjects(y'); % subjects
    x = iC;
    % standard error of mean
    err = std(y)/sqrt(numel(y));

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
                errorbar(x,mean(y,'omitnan'),err,'Marker','.','MarkerSize',style.scatter.MarkerSize,'MarkerFaceColor',faceColor,'MarkerEdgeColor',faceColor,...
        'Color',faceColor,'LineWidth',2);
    end

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
ylabel(sprintf('Slope (\\Delta ITPC/s)'))

xlabel('Precue') % need to retain sizing?
xticks([1])
xticklabels({'All'})
xlim([1-style.xBuffer/1.5 1+style.xBuffer/1.5])
ylim([0 0.1])

% Stats annotation
if annotateStats
    txt = meg_annotateStats(1,max(fh.YLim)*0.97,'**'); 
end

if saveFigs
    figTitle = sprintf('meg_manuscriptFigs_TE_allTrials_mdlFit_slope_bar_%s',dateStr);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end

%% Export csv long format data for R 
exportCSV = 1;
cueLevel = {'all'}; 
clear V T 
tableHeaders = {'subject','session','precue','dummy'...
    'intercept','slope'};
varNames = {'intercept','slope'}; 
count = 1;
for iS = 1:10 % subjects
    for iSession = 1:2
        for iC = 1:numel(cueLevel)
            for iD = 1:2
                for iV = 1:numel(varNames)
                    V.subject(count) = iS;
                    V.session(count) = iSession;
                    V.precue(count) = iC;
                    V.dummy(count) = iD-1; 
                    if iSession==1
                        idxS = iS*2-1;
                    elseif iSession==2
                        idxS = iS*2;
                    end
                    if iD==1 % if dummy 
                        V.(tableHeaders{iV+4})(count) = 0; % slope to 0 
                    else 
                        V.(tableHeaders{iV+4})(count) = mdlFit.linear.(cueLevel{iC}).session.minSolution(idxS,iV);
                    end
                end
                count = count+1;
            end
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
    csvName = sprintf('TANoise_MdlFit_Linear_allTrials');
    csvPath = sprintf('%s/%s.csv', csvDir, csvName);
    writetable(T,csvPath)
end

%% ttest
% Average to subjects
for iV = 1:2
    mdlFit.linear.all.session.minSolutionSubjects(:,iV) = meg_sessions2subjects(mdlFit.linear.all.session.minSolution(:,iV)'); 
end

[test.h test.p test.ci test.stats] = ttest(mdlFit.linear.all.session.minSolutionSubjects(:,2));


