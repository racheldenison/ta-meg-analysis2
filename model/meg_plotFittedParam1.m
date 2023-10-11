function meg_plotFittedParam1(fittedP,p,fitType,figDir)
% structure: each field corresponds to one fitted parameter 
% 10 sub x 2 sessions x 2 precues

%% Settings
plotSubjects = 1; 
xJitter = 0.2; 
idxPeaks = 1:10; 

%% Setup 
fields = fieldnames(fittedP); 
nFields = numel(fields);

xs = repmat(1,[10,1]);
sSize = 100; 
subjectColor = [0.7 0.7 0.7]; 

%% Saving
saveFigs = 1; % 0 
figFormat = 'png'; 

%% --- Figure (sessions) ---
for iV = 1:nFields
    figure
    set(gcf,'Position',[100 100 300 300])
    hold on
    
    % -- All other params ---
    for iS = 1:2 % sessions
        if iS == 1
            scatter(fittedP.(fields{iV})(:,iS,1),... % precue T1 
                fittedP.(fields{iV})(:,iS,2),... % precue T2 
                sSize, 'filled','MarkerFaceColor','k','MarkerEdgeColor','w')
        else
            scatter(fittedP.(fields{iV})(:,iS,1),...
                fittedP.(fields{iV})(:,iS,2),...
                sSize, 'filled','MarkerFaceColor','w','MarkerEdgeColor','k')
        end
    end

    for iSubject = 1:10
        % --- Subject lines (between S1 and S2) ---
        s = plot([fittedP.(fields{iV})(iSubject,1,1) fittedP.(fields{iV})(iSubject,2,1)],...
            [fittedP.(fields{iV})(iSubject,1,2) fittedP.(fields{iV})(iSubject,2,2)],...
            'Color',subjectColor);
    end

    % --- Ref line --- 
    r = refline(1,0);
    r.Color = 'k'; 
    r.LineStyle = ':';
    r.LineWidth = 1; 

    % --- Format ---
    title(und2space(sprintf('%s',fields{iV})))
    xlabel('Precue T1')
    ylabel('Precue T2')
    meg_figureStyle
    axis square 
    axis equal 

    % --- Save fig ---
    if saveFigs
        figTitle = sprintf('TANoise_ITPCFit_FittedParams_sessions_%s_%s',fitType,und2space(fields{iV}));
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle,figFormat))
    end
end

%% --- Figure (subjects) ---
for iV = 1:nFields
    figure
    set(gcf,'Position',[100 100 300 300])
    hold on
    
    % -- All other params ---
    scatter(mean(fittedP.(fields{iV})(:,:,1),2),... % precue T1
        mean(fittedP.(fields{iV})(:,iS,2),2),... % precue T2
        sSize, 'filled','MarkerFaceColor',subjectColor,'MarkerEdgeColor','w')

    % --- Ref line --- 
    r = refline(1,0);
    r.Color = 'k'; 
    r.LineStyle = ':';
    r.LineWidth = 1; 

    % --- Format ---
    title(und2space(sprintf('%s',fields{iV})))
    xlabel('Precue T1')
    ylabel('Precue T2')
    meg_figureStyle
    axis square 
    axis equal 

    % --- Save fig ---
    if saveFigs
        figTitle = sprintf('TANoise_ITPCFit_FittedParams_subjects_%s_%s',fitType,und2space(fields{iV}));
        saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle,figFormat))
    end
end

%% Polar plot of phase (sessions) 
for iV = 1:nFields
    if strcmp(fields{iV},'phaserad')
                figure
        set(gcf,'Position',[100 100 300 300])
        for iSession = 1:2
            theta = circ_dist(fittedP.(fields{iV})(:,iSession,1),fittedP.(fields{iV})(:,iSession,2));
            rho = repmat(1,[10,1]);
            if iSession == 1 
                polarscatter(theta,rho,sSize,'filled','MarkerFaceColor','k','MarkerEdgeColor','w')
            else
                polarscatter(theta,rho,sSize,'filled','MarkerFaceColor','w','MarkerEdgeColor','k')
            end
            hold on 
        end

        % --- Color by mean amplitude across precue conditions ---
        for iS = 1:10
            for iSession = 1:2
                theta = circ_dist(fittedP.(fields{iV})(iS,iSession,1),fittedP.(fields{iV})(iS,iSession,2));
                rho = 2;
                amp = mean(fittedP.amplitude(iS,iSession,:),3); 
                
                edges = 0:0.01:0.05; 
                bin = discretize(amp,edges); 

                ampMap = [204 206 251;
                    164 167 227;
                    127 127 203;
                    89 90 180;
                    50 52 157;
                    10 13 134]/255;

                scat = polarscatter(theta,rho,sSize,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
                scat.MarkerFaceColor = ampMap(bin,:);
            end
        end
        
    end
end
% --- Save fig ---
if saveFigs
    figTitle = sprintf('TANoise_ITPCFit_FittedParams_sessions_%s_phaseRad',fitType);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle,figFormat))
end

%% Polar plot of phase (subjects) 
for iV = 1:nFields
    if strcmp(fields{iV},'phaserad')
        figure
        set(gcf,'Position',[100 100 300 300]) 
            theta = circ_dist(circ_mean(fittedP.(fields{iV})(:,:,1),[],2),...
                circ_mean(fittedP.(fields{iV})(:,:,2),[],2));
            rho = repmat(1,[10,1]);
            polarscatter(theta,rho,sSize,'filled','MarkerFaceColor',subjectColor)
            hold on 

            % --- Color by mean amplitude across precue conditions ---
            % for iS = 1:10
            %     for iSession = 1:2
            %         theta = circ_dist(fittedP.(fields{iV})(iS,iSession,1),fittedP.(fields{iV})(iS,iSession,2));
            %         rho = 2;
            %         amp = mean(fittedP.amplitude(iS,iSession,:),3);
            % 
            %         edges = 0:0.01:0.05;
            %         bin = discretize(amp,edges);
            % 
            %         ampMap = [204 206 251;
            %             164 167 227;
            %             127 127 203;
            %             89 90 180;
            %             50 52 157;
            %             10 13 134]/255;
            % 
            %         scat = polarscatter(theta,rho,sSize,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
            %         scat.MarkerFaceColor = ampMap(bin,:);
            %     end
            % end


    end
end
% --- Save fig ---
if saveFigs
    figTitle = sprintf('TANoise_ITPCFit_FittedParams_subjects_%s_phaseRad',fitType);
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle,figFormat))
end

%% Plot data by precue 
% figure
% set(gcf,'Position',[100 100 1200 300])
% for iV = 1:nFields
%     subplot(1,nFields,iV)
%     hold on
% 
%     % -- All other params ---
%     for iC = 1:2 % precues
%         for iS = 1:2 % sessions
%             if iS == 1
%                 sign = -1;
%             else
%                 sign = 1;
%             end
%             scatter( (sign*xJitter)+xs*iC,...
%                 fittedP.(fields{iV})(:,iS,iC),...
%                 sSize, 'filled','MarkerFaceColor',p.cueColors(iC,:))
%         end
% 
%         for iSubject = 1:10
%             % --- Subject lines (between S1 and S2) ---
%             s = plot([(-1*xJitter)+xs(1)*iC (1*xJitter)+xs(1)*iC],...
%                 [fittedP.(fields{iV})(iSubject,1,iC) fittedP.(fields{iV})(iSubject,2,iC)],...
%                 'Color',subjectColor);
%         end
%     end
% 
%     % --- Format --- 
%     ylabel(sprintf('%s',fields{iV}))
%     meg_figureStyle
% end





