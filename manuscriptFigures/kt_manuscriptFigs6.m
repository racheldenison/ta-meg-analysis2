%% FIG: 20 Hz ITPC peaks T1 and T2 by precue 
% December 20, 2021 
% something isn't right w these vals... 

%% Load data and setup 
load('/Users/kantian/Dropbox/Data/TANoise/fromRachel/itpcNorm_TS_Peaks_N10_20211225_workspace.mat')
p = meg_params('TANoise_ITPCsession8'); 
plotSubjects = 1; 

%% Fig setup 
% Normalized ITPC 
% att (precue T1, precue T2 x target (T1, T2) x subject (10) x session (2) 
peakDataAve = mean(peakData,4); % average sessions 

xJitter = 0.2; 
idxPeaks = 1:10; % 10 subjects 
idxPeaks(5) = []; % select subjects (9) 

sSize = 50; % subject dots marker size 
subjectColor = [0.86 0.86 0.86]; % subject lines 

peak_T1 = squeeze(peakDataAve(:,1,:)); 
peak_T2 = squeeze(peakDataAve(:,2,:));

%% Plot normalized ITPC (subject dots) 
% --- Figure ---
figure
set(gcf,'Position',[100 100 300 400])
hold on 
meg_figureStyle

% T1 
% PrecueT1 
x_valid = repmat(1-xJitter,[9,1]);
% standard error of difference 
errorbar(x_valid(1),mean(peak_T1(1,idxPeaks)),peakDiffSte(1),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(1,:),'MarkerEdgeColor',p.cueColors(1,:),...
    'Color',p.cueColors(1,:),'LineWidth',2);
if plotSubjects 
    scatter(x_valid,peak_T1(1,idxPeaks),'filled','MarkerFaceColor','w')
    scatter(x_valid,peak_T1(1,idxPeaks),sSize,'filled','MarkerFaceColor',p.cueColors(1,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
end
% PrecueT2 
x_invalid = repmat(1+xJitter,[9,1]);
errorbar(x_invalid(1),mean(peak_T1(2,idxPeaks)),peakDiffSte(1),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(2,:),'MarkerEdgeColor',p.cueColors(2,:),...
    'Color',p.cueColors(2,:),'LineWidth',2);
if plotSubjects
    % subject lines 
    for i = idxPeaks
        s = plot([x_valid x_invalid],[peak_T1(1,i) peak_T1(2,i)],'Color',subjectColor);
    end
    scatter(x_invalid,peak_T1(2,idxPeaks),'filled','MarkerFaceColor','w')
    scatter(x_invalid,peak_T1(2,idxPeaks),sSize,'filled','MarkerFaceColor',p.cueColors(2,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
end

% T2 
% PrecueT1 
x_valid = repmat(2-xJitter,[9,1]);
errorbar(x_valid(1),mean(peak_T2(1,idxPeaks)),peakDiffSte(2),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(1,:),'MarkerEdgeColor',p.cueColors(1,:),...
    'Color',p.cueColors(1,:),'LineWidth',2);
if plotSubjects 
    scatter(x_valid,peak_T2(1,idxPeaks),'filled','MarkerFaceColor','w')
    scatter(x_valid,peak_T2(1,idxPeaks),'filled','MarkerFaceColor',p.cueColors(1,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
end
% PrecueT2 
x_invalid = repmat(2+xJitter,[9,1]);
errorbar(x_invalid(1),mean(peak_T2(2,idxPeaks)),peakDiffSte(2),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(2,:),'MarkerEdgeColor',p.cueColors(2,:),...
    'Color',p.cueColors(2,:),'LineWidth',2);
if plotSubjects 
    % subject lines 
    for i = idxPeaks
        s = plot([x_valid x_invalid],[peak_T2(1,i) peak_T2(2,i)],'Color',subjectColor);
    end
    scatter(x_invalid,peak_T2(2,idxPeaks),'filled','MarkerFaceColor','w')
    scatter(x_invalid,peak_T2(2,idxPeaks),'filled','MarkerFaceColor',p.cueColors(2,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
end

xlim([0 3])
xticks([1,2])
xticklabels({'T1','T2'}); 
ylabel('Normalized ITPC')

% lines between T1 and T2 
if plotSubjects
    % subject lines 
    for i = idxPeaks
        s = plot([repmat(1+xJitter,[9,1]) repmat(2-xJitter,[9,1])],[peak_T1(2,i) peak_T2(1,i)],'Color',subjectColor);
    end
end


