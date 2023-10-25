% plot TA2 ERF and peaks by top 5 visually responsive channels

%% Settings
plotFigs = 1; 
saveFigs = 1; 
plotNeutral = 0; 

%% Load data
addpath(genpath('../../ta-meg-analysis2/'))

% load data D structure
load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/figures/cue_max/nCh_1_5/ERFbyCue.mat')
% A structure
channelGroups = [1,2,3,5,7,10,15,20,50];
for i = 1:numel(channelGroups)
    vals = load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/Group/figures/cue_max/nCh_1_%d/max.mat',channelGroups(i)));
    ACh(i) = vals;
end

% Get parameters
p = meg_params('TANoise_Analysis');

%% Figure style
color_pT1 = [142 195 221]/255; % blue
color_pT2 = [189 76 120]/255; % pink 
color_pN  = [126 126 126]/255;
color_grey = [0.6 0.6 0.6];

lineWidth = 2;
unit = 10e-15; % femto unit
alpha = 0.8;
markerSize = 20;
x_attOffset = 0.3; 

%% Plot ERF
if plotFigs
    figure
    set(gcf,'Position',[0 0 900 800])

    h(1) = subplot (3,2,[1 2]); 
    h(1).Position(3) = h(1).Position(3) * 0.7; 
    meg_figureStyle
    hold on
    t = p.tstart:p.tstop;
    tPlot = find(t==-50):10:find(t==2300);

    for iK = 1:numel(p.eventTimes)
        xline(p.eventTimes(iK),'Color',color_grey)
    end

    yline(0,'Color',color_grey)

    y_precueT1 = D.group.cueT1(tPlot)/unit;
    y_precueT2 = D.group.cueT2(tPlot)/unit;
    if plotNeutral
        y_precueN = D.group.cueN(tPlot)/unit;
        pn = plot(t(tPlot),y_precueN,'LineWidth',lineWidth,'color',color_pN);
        pn.Color(4) = alpha;
    end

    p1 = plot(t(tPlot),y_precueT1,'LineWidth',lineWidth,'color',color_pT1);
    p2 = plot(t(tPlot),y_precueT2,'LineWidth',lineWidth,'color',color_pT2);

    p1.Color(4) = alpha;
    p2.Color(4) = alpha;

    xlims = [-50,2300];
    ylims = [-8,18];
    xlim(xlims)
    ylim(ylims)
    line([xlims(1) xlims(2)],[ylims(2) ylims(2)],'Color','k')
    line([xlims(2) xlims(2)],[ylims(1) ylims(2)],'Color','k')

    xlabel('Time (ms)')
    ylabel('Amplitude (fT)')

    legend([p1,p2],'Precue T1','Precue T2')
    legend boxoff  

    % Plot peak latency
    targets = {'T1','T2'};
    if plotNeutral
       precues = {'cueT1','cueT2','cueN'};
    else
       precues = {'cueT1','cueT2'}; 
    end

    for iT = 1:numel(targets)
        subplot (3,2,2+iT)
        meg_figureStyle
        hold on
        ylabel(sprintf('%s peak\n latency (ms)',targets{iT}))
        ylim([100 150])
        xlim([0 numel(channelGroups)+1])
        xticks(1:numel(channelGroups))
        xticklabels(channelGroups)
        if plotNeutral 
            nAtt = 3; 
        else
            nAtt = 2; 
        end
        for iCh = 1:numel(channelGroups)
            for iP = 1:nAtt
                n = numel(ACh(iCh).A.(targets{iT}).(precues{iP}).maxIdx);
                x = (iCh - x_attOffset*2) + (x_attOffset * iP);
                y = ACh(iCh).A.(targets{iT}).(precues{iP}).maxIdx - p.eventTimes(iT+1);
                y = mean(y,'omitnan');
                err = ACh(iCh).A.(targets{iT}).(precues{iP}).maxIdxStd/sqrt(n);
                e = errorbar(x, y, err,'Marker','.','MarkerSize',markerSize);
                e.CapSize = 0;
                e.LineWidth = 2;
                if iP==1
                    e.Color = color_pT1;
                elseif iP==2
                    e.Color = color_pT2;
                elseif iP==3
                    e.Color = color_pN;
                end
            end
        end
    end

    % Plot peak amplitude
    for iT = 1:2
        subplot (3,2,4+iT)
        meg_figureStyle
        hold on
        ylabel(sprintf('%s peak\n amplitude (fT)',targets{iT}))
        xlim([0 numel(channelGroups)+1])
        xticks(1:numel(channelGroups))
        xticklabels(channelGroups)
        ylim([0 30])
        if iT == 2
            xlabel('Number of channels')
        end
        for iCh = 1:numel(channelGroups)
            for iP = 1:nAtt
                n = numel(ACh(iCh).A.(targets{iT}).(precues{iP}).maxIdx);
                x = (iCh - x_attOffset*2) + (x_attOffset * iP);
                y = ACh(iCh).A.(targets{iT}).(precues{iP}).max/unit;
                y = mean(y,'omitnan');
                err = (ACh(iCh).A.(targets{iT}).(precues{iP}).maxStd/unit)/sqrt(n);
                e = errorbar(x, y, err,'Marker','.','MarkerSize',markerSize);
                e.CapSize = 0;
                e.LineWidth = 2;
                if iP==1
                    e.Color = color_pT1;
                elseif iP==2
                    e.Color = color_pT2;
                elseif iP==3
                    e.Color = color_pN;
                end
            end
        end
    end

    if saveFigs
        % Save
        figName = 'TA2_Peak_ERF';
        figFolder = '/Users/kantian/Dropbox/Data/TA2/ManuscriptFigures';
        saveas(gcf,sprintf('%s/%s.svg', figFolder, figName))
    end
end

%% Save 
save('TA2_MaxLatency_TargetPrecue.mat','ACh','channelGroups','D','-v7.3'); 

%% Compile top 20 sessions
for iSession = 1:20
    channels(iSession,:) = groupC(iSession).sortChByProm; 
end

save('TA2_channels_rankPkProm.mat','channels','-v7.3')

