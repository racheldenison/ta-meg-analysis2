function meg_manuscriptFigs_TE_allTrials_mdlFit
% Plot ITPC all trials with model fits (linear, and linear+periodic)

%% Load data
load('/Users/kantian/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ModelFit_1_freeFreq_GroupTS_231102.mat')

%% Figure settings 
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotErrorBars = 1; % turn off for subject-level plots 
restrictYLim = 0; % turn on for group-level manuscript matching ylims 
plotPrePrecue = 0; % bar showing pre-precue baseline period 
saveFigs = 1; 
cueLevel = {'all'}; 

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
        [style, colors] = meg_manuscriptStyle;
        xlim([-300 2400]); % xlim([-100 2400])
        if restrictYLim
            ylim([0.25 0.4])
            yticks(0.25:0.05:0.4)
        end
        xlabel('Time (ms)')
        ylabel('ITPC')

        for iC = 1:numel(cueLevel)
            % --- Plot data error bars ---
            % Requires https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
            sampling = 1:10:7001;
            errorBarType = 'EB-TS-baselineNorm'; % 'EB-TS' 'EB-fit' 'EB-TS-baselineNorm'
            if plotErrorBars
                includeIdx = 1:10;
                groupVals = data.(cueLevel{iC}).group;
                vals = data.(cueLevel{iC}).subject;
                switch errorBarType
                    % case 'EB-TS'
                    %     t = p.t(sampling);
                    %     line_allTrials = shadedErrorBar(t, groupVals(sampling), std(vals(sampling,includeIdx),[],2)/sqrt(numel(includeIdx)),...
                    %         'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',1);
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
                            line_allTrials.edge(iL).LineWidth = 0.2;
                        end
                        % case 'EB-fit'
                        %     t = p.t(toi);
                        %     % Collect subject-level fits into single matrix
                        %     for i = 1:numel(c_subs)
                        %         sub_fits(:,i) = c_subs(i).y1_est;
                        %     end
                        %     line_allTrials = shadedErrorBar(t, y1_est, std(sub_fits(:,includeIdx),[],2)/sqrt(numel(includeIdx)),...
                        %         'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',1);
                end
            end

            % --- Data (by precue) ---
            dataFit = data.(cueLevel{iC}).(fitLevel)(:,iS)';

            % --- Plot data ---
            plot(p.t(sampling),dataFit(sampling),'LineWidth',2,'Color',colors.mediumPurple)

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
                    paramNames = {'intercept','slope','amplitude','phase','freq'};
                    % --- Annotate fitted frequency ---
                    fStr = sprintf('%0.2f Hz',fittedX(:,:,strcmp(paramNames,'freq')) );
                    fStrTxt = text(30+p.eventTimes(1),yl(1)+0.003,fStr,'HorizontalAlignment','left','VerticalAlignment','bottom');
                    fStrTxt.FontSize = 14;
                    fStrTxt.FontName = 'Helvetica-Light';

                case 'linear'
                    paramNames = {'intercept','slope'};
            end
            [~,yhat] = meg_objectiveFunction1_freeFreq(fittedX,dataFit(toi),t,Fs,paramNames,fitTypes{iF});
            plot(t,yhat,'--','LineWidth',2,'Color',colors.darkPurple)
        end

        % --- Plot shaded baseline window ---
        xl = xlim;
        yl = ylim;
        x = [p.t(tIdx(1)) p.t(tIdx(end)) p.t(tIdx(end)) p.t(tIdx(1))];
        % y = [yl(1) yl(1) yl(2) yl(2)]; % [yl(1) yl(1) yl(2) yl(2)]
        % patch(x,y,lightgrey,'EdgeColor',lightgrey)
        yScale = yl(2)-yl(1);
        ySize = (0.02*0.15)*(yScale/0.15);
        y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize]; % [yl(1) yl(1) yl(2) yl(2)]

        patch(x,y,colors.mediumgrey,'EdgeColor',colors.mediumgrey)
        % plot([x(1) x(2)],[y(1) y(1)],'Color',mediumgrey,'LineWidth',30)

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

        % --- Plot event lines ---
        p.eventNamesCap = {'Precue','T1','T2','Response cue'};
        for i = 1:numel(p.eventTimes)
            xh = xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1);
            meg_sendToBack(xh)
            hold on
            % --- Plot event names ---
            if i==1 % Color the precue
                % Precue (all trials)
                ySet = max(fh.YLim)+(diff(fh.YLim)*0.01); 
                txtName = sprintf('Precue\n(All trials)'); 
                txt = text(p.eventTimes(i),ySet,txtName,'EdgeColor','none',...
                    'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Bottom');
                txt.Color = colors.darkPurple;
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
            figTitle = sprintf('meg_manuscriptFigs_TE_allTrials_mdlFit_%s_%s',fitTypes{iF},dateStr); 
            saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
        end
    end
end