function meg_manuscriptFigs_TE_byPrecue_mdlFit
% function meg_manuscriptFigs_TE_byPrecue_mdlFit
% Plot 2 Hz model fit 

%% Load data
user = 'kantian'; 
% Model fit data 
filename = sprintf('/Users/kantian/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ModelFit_SeparatePrecueT1T2_2Hz_GroupTS_231218.mat',user);
% ITPC data 
load(filename)

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotErrorBars = 1; % turn off for subject-level plots 
restrictYLim = 1; % turn on for group-level manuscript matching ylims 
plotPrePrecue = 0; % bar showing pre-precue baseline period 
saveFigs = 1; 
cueLevel = {'cueT1','cueT2'}; 
[style, colors] = meg_manuscriptStyle;

% Figure directory 
user = 'kantian'; 
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/manuscriptFigures/figs',user); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Plot figure
for iF = 1:numel(fitTypes)
    for iS = 1:size(data.all.(fitLevel),2) % sessions
        figure
        switch fitLevel
            case {'session','subject','group'}
                set(gcf,'Position',[100 100 500 300]) %  [100 100 600 400]
            otherwise
                error('Specify session- or subject-level fit')
        end
        fh = subplot(1,1,1);
        hold on

        % --- Format ---
        meg_figureStyle
        xlim([-300 2400]); % xlim([-100 2400])
        if restrictYLim
            ylim([0.26 0.4])
            yticks(0.26:0.02:0.4)
        end
        xlabel('Time (ms)')
        ylabel('ITPC')

        for iC = 1:numel(cueLevel)
            % --- Set cue colors ---
            if iC==2
                cLight = 'lightRed';
                cMed = 'mediumRed';
            elseif iC==1
                cLight = 'lightBlue';
                cMed = 'mediumBlue';
            end

            % --- Plot data error bars ---
            % Requires https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
            sampling = 1:10:7001;
            errorBarType = 'SED'; % 'EB-TS' 'EB-fit' 'EB-TS-baselineNorm'
            if plotErrorBars
                includeIdx = 1:10;
                groupVals = data.(cueLevel{iC}).group;
                vals = data.(cueLevel{iC}).subject;
                switch errorBarType
                    case 'SED'
                        ydiff = data.(cueLevel{1}).subject(:,includeIdx)-data.(cueLevel{2}).subject(:,includeIdx);
                        sed = std(ydiff,[],2,'omitnan')./sqrt(numel(includeIdx));
                        % Plot
                        t_plot = p.t(sampling);
                        line_allTrials = shadedErrorBar(t_plot, groupVals(sampling), sed(sampling),...
                            'lineProps', {'MarkerFaceColor',p.cueColors(iC,:), 'LineWidth', 0.2, 'Color',p.cueColors(iC,:)}, 'transparent',1);
                        % line_allTrials.patch.FaceColor = colors.(cLight); p.cueColors(iC,:)
                        line_allTrials.patch.FaceAlpha = 0.3;
                        for iL = 1:2
                            line_allTrials.edge(iL).Color = colors.(cLight);
                            line_allTrials.edge(iL).LineWidth = style.ebLineWidth;
                        end
                    case 'EB-TS-baselineNorm'
                        baselineTOI = -100:-1; % why is this -200?
                        baselineTOIIdx = find(p.t==baselineTOI(1)):find(p.t==baselineTOI(end));
                        for i = 1:size(vals,2) % average per subject
                            baselineS(i) = mean(vals(baselineTOIIdx,i),1);
                            normVals(:,i) = vals(:,i)./baselineS(i);
                        end
                        % sem in units of percent change
                        sem = std(normVals(sampling,includeIdx),[],2)/sqrt(numel(includeIdx));
                        % scale sem back to ITPC units
                        baselineG = mean(groupVals(baselineTOIIdx));
                        semScaled = sem*baselineG;
                        % Plot
                        t_plot = p.t(sampling);
                        line_allTrials = shadedErrorBar(t_plot, groupVals(sampling), semScaled,...
                            'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',0);
                        line_allTrials.patch.FaceColor = colors.lightPurple;
                        line_allTrials.patch.FaceAlpha = 1;
                        for iL = 1:2
                            line_allTrials.edge(iL).Color = colors.mediumPurple;
                            line_allTrials.edge(iL).LineWidth = style.ebLineWidth;
                        end
                end
            end
        end

        % --- Plot shaded baseline window ---
        xl = xlim;
        yl = ylim;
        x = [p.t(tIdx(1)) p.t(tIdx(end)) p.t(tIdx(end)) p.t(tIdx(1))];
        yScale = yl(2)-yl(1);
        ySize = (0.02*0.15)*(yScale/0.15);
        y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize]; % [yl(1) yl(1) yl(2) yl(2)]

        patch(x,y,colors.mediumgrey,'EdgeColor',colors.mediumgrey)

        % --- Plot pre-precue baseline check ---
        if plotPrePrecue
            x = [p.t(tIdx(1))-100 p.t(tIdx(1)) p.t(tIdx(1)) p.t(tIdx(1))-100];
            y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize];
            patch(x,y,colors.lightgrey,'EdgeColor',colors.lightgrey)
        end

        if showN
            % --- Add n annotation ---
            nStr = sprintf('n = %d x 2 sessions',size(data.all.subject,2));
            nStrTxt = text(0.945*xl(2),yl(1)+0.003,nStr,'HorizontalAlignment','right','VerticalAlignment','bottom');
            nStrTxt.FontSize = 14;
            nStrTxt.FontName = 'Helvetica-Light';
        end

        % set(gca,'children',flipud(get(gca,'children'))) % change order of plotted objects
        set(gca ,'Layer', 'Top') % bring axes to front

        for iC = 1:numel(cueLevel) % Plot again to be on top
            % --- Set cue colors ---
            if iC==2
                cLight = 'lightRed';
                cMed = 'mediumRed';
            elseif iC==1
                cLight = 'lightBlue';
                cMed = 'mediumBlue';
            end

            % --- Data (by precue) ---
            dataFit = data.(cueLevel{iC}).(fitLevel)(:,iS)';

            % --- Plot data ---
            plot(p.t(sampling),dataFit(sampling),'LineWidth',2,'Color',p.cueColors(iC,:))

            % --- Plot model fit ---
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS)  );
            if iC == 1
                idx1 = idx;
            elseif iC ==2
                idx2 = idx;
            end
            fittedX = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:);
            switch fitTypes{iF}
                case 'linear2Hz'
                    paramNames = {'intercept','slope','amplitude','phase'};
                case 'linear'
                    paramNames = {'intercept','slope'};
            end
            [~,yhat] = meg_objectiveFunction1(fittedX,dataFit(toi),t,Fs,paramNames,fitTypes{iF});
            plot(t,yhat,'--','LineWidth',2,'Color',colors.(cMed))
        end

        % --- Annotate fitted frequency ---
        switch fitTypes{iF}
            case 'linear2Hz'
                fStr = sprintf('%d Hz',2);
                fStrTxt = text(30+p.eventTimes(1),yl(1)+0.003,fStr,'HorizontalAlignment','left','VerticalAlignment','bottom');
                fStrTxt.FontSize = 14;
                fStrTxt.FontName = 'Helvetica-Light';
                fStrTxt.Color = 'k';
        end

        % --- Plot event lines ---
        p.eventNamesCap = {'Precue','T1','T2','Response cue'};
        % p.eventNamesCap{1} = sprintf('Precue T1\nPrecue T2');
        for i = 1:numel(p.eventTimes)
            xh = xline(p.eventTimes(i),'Color',colors.eventLines,'LineWidth',1);
            meg_sendToBack(xh)
            hold on
            % --- Plot event names ---
            if i==1 % Color the precue
                % Precue T1
                yOffset = diff(fh.YLim)*0.07;
                ySet = max(fh.YLim)+(diff(fh.YLim)*0.01);
                txt = text(p.eventTimes(i),ySet+yOffset,'Precue T1','EdgeColor','none',...
                    'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Bottom');
                txt.Color = p.cueColors(1,:); % colors.mediumBlue;
                % Precue T2
                txt = text(p.eventTimes(i),ySet,'Precue T2','EdgeColor','none',...
                    'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Bottom');
                txt.Color = p.cueColors(2,:); % colors.mediumRed;
            elseif i==4 % Response cue
                text(p.eventTimes(i),ySet,p.eventNamesCap{i},'EdgeColor','none',...
                    'FontSize',14,'HorizontalAlignment','right','VerticalAlignment','Bottom');
            else
                text(p.eventTimes(i),ySet,p.eventNamesCap{i},'EdgeColor','none',...
                    'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','Bottom');
            end
        end

        % --- Titles --
        if titleVis
            switch fitTypes{iF}
                case 'linear2Hz'
                    titleText_Fitted = sprintf('Fitted: intercept = %0.2f, slope = %0.3f, amplitude = %0.2f, phase = %0.2f\\pi\nfreq = %0.2f',...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2),...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,3),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,4)/pi,...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,5) );
                case 'linear'
                    titleText_Fitted = sprintf('Fitted: intercept = %0.2f, slope = %0.3f',...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2));
            end
            switch fitLevel
                case {'sessions','subjects'}
                    titleText1 = sprintf('Fits to all trials (%s), fval = %0.2e, %s',...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(iS),...
                        und2space(sessionNames{iS}));
                case 'group'
                    titleText1 = sprintf('Fits to all trials (%s), fval = %0.2e, n = %d',...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
                        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(iS),...
                        size(data.all.subject,2));
            end
            titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted);
            title(titleText)

            ax = gca;
            ax.TitleFontSizeMultiplier = 0.6;
            ax.TitleFontWeight = 'normal';
        end

        % --- Save fig ---
        if saveFigs
            figTitle = sprintf('meg_manuscriptFigs_TE_byPrecue_mdlFit_%s_%s',fitTypes{iF},dateStr); 
            saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
        end
    end
end




