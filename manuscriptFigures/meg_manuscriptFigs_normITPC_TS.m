% Plot normalized ITPC TS by att w shaded error bars 

%% Settings
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotErrorBars = 1; % turn off for subject-level plots 
restrictYLim = 0; % turn on for group-level manuscript matching ylims 
saveFigs = 1; 

% Figure directory 
user = 'kantian'; 
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/manuscriptFigures/figs',user); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% MEG settings 
p = meg_params('TANoise_ITPCsession8'); 
cueLevel = {'cueT1','cueT2'};
includeIdx = [1,2,3,4,6,7,8,9,10];

% --- Timing ---
sampling = 1:10:7001; % softens curves
foi = 20; % frequency of interest, Hz
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 

%% Load data 
filename = sprintf('/Users/%s/Dropbox/Data/TANoise/MEG/Group/mat/groupA_ITPCspectrogram_byAtt.mat',user); 
load(filename)

%% Plot normalized ITPC time series with shaded error bars 
figure;
set(gcf,'Position',[100 100 500 300]) %  [100 100 600 400]
fh = subplot(1,1,1);
hold on

% --- Format ---
meg_figureStyle
[style, colors] = meg_manuscriptStyle;
xlim([-300 2400]); % xlim([-100 2400])
xlabel('Time (ms)')
ylabel('Normalized ITPC')

% --- Plot 0 line ---
xh = yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1); 
meg_sendToBack(xh)

% --- Plot data error bars ---
% Requires https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
errorBarType = 'SED';
for iC = 1:numel(cueLevel)
    % --- Get data ---
    groupVals = mean(A.(cueLevel{iC}).normSubjectFlipped(20,:,includeIdx),3,'omitnan');
    % --- Set cue colors ---
    if iC==2
        cLight = 'lightRed';
        cMed = 'mediumRed';
    elseif iC==1
        cLight = 'lightBlue';
        cMed = 'mediumBlue';
    end

    if plotErrorBars
        switch errorBarType
            case 'SED' % standard error of the difference
                ydiff = squeeze(A.cueT1.normSubjectFlipped(20,:,includeIdx))-squeeze(A.cueT2.normSubjectFlipped(20,:,includeIdx));
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
        end
    end
end

for iC = 1:numel(cueLevel)
    groupVals = mean(A.(cueLevel{iC}).normSubjectFlipped(20,:,includeIdx),3,'omitnan');
    % --- Plot data (on top) ---
    plot(p.t(sampling),groupVals(sampling),'LineWidth',2,'Color',p.cueColors(iC,:))
end

% --- Plot shaded baseline window ---
xl = xlim;
yl = ylim;
x = [p.t(tIdx(1)) p.t(tIdx(end)) p.t(tIdx(end)) p.t(tIdx(1))];
yScale = yl(2)-yl(1);
ySize = (0.02*0.15)*(yScale/0.15);
y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize]; % [yl(1) yl(1) yl(2) yl(2)]

patch(x,y,colors.mediumgrey,'EdgeColor',colors.mediumgrey)

if showN
    % --- Add n annotation ---
    nStr = sprintf('n = %d x 2 sessions',size(includeIdx,2));
    nStrTxt = text(0.945*xl(2),yl(1)+0.003,nStr,'HorizontalAlignment','right','VerticalAlignment','bottom');
    nStrTxt.FontSize = 14;
    nStrTxt.FontName = 'Helvetica-Light';
end

% --- Plot event lines ---
p.eventNamesCap = {'Precue','T1','T2','Response cue'};
% p.eventNamesCap{1} = sprintf('Precue T1\nPrecue T2');
for i = 1:numel(p.eventTimes)
    xh = xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1);
    meg_sendToBack(xh)
    hold on
    % --- Plot event names ---
    if i==1 % Color the precue
        % Precue T1
        yOffset = diff(fh.YLim)*0.07;
        ySet = max(fh.YLim)+(diff(fh.YLim)*0.01); 
        txt = text(p.eventTimes(i),ySet+yOffset,'Precue T1','EdgeColor','none',...
            'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Bottom');
        txt.Color = colors.mediumBlue;
        % Precue T2
        txt = text(p.eventTimes(i),ySet,'Precue T2','EdgeColor','none',...
            'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','Bottom');
        txt.Color = colors.mediumRed;
    elseif i==4 % Response cue
        text(p.eventTimes(i),ySet,p.eventNamesCap{i},'EdgeColor','none',...
            'FontSize',14,'HorizontalAlignment','right','VerticalAlignment','Bottom');
    else
        text(p.eventTimes(i),ySet,p.eventNamesCap{i},'EdgeColor','none',...
            'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','Bottom');
    end
end

% --- Save fig ---
if saveFigs
    figTitle = sprintf('meg_manuscriptFigs_normITPC_TS_%s',dateStr);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
end






