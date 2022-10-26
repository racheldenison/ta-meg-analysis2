function meg_fig7

% a. Plots group normalized ITPC 
% b. Plots normalized 20 Hz ITPC peak target-evoked response for each target and precue condition 

%% Load data and setup 
variable = 'Average trial power'; % ITPC
pks = load('ITPC/peaks.mat'); % load peak variables
pks.pkVals = []; 
pks.pkVals = pkVals; 
pks.pkHeaders = {'subject','session','cue','target'}; 

p = meg_params('TANoise_ITPCsession8'); 
plotSubjects = 1; 

%% Fig setup 
% Normalized ITPC 
% att (precue T1, precue T2 x target (T1, T2) x subject (10) x session (2) 
peakDataAve = mean(pks.pkVals,4); % average sessions 

xJitter = 0.2; 
idxPeaks = 1:10; % 10 subjects 
idxPeaks(5) = []; % select subjects (9) 

sSize = 50; % subject dots marker size 
subjectColor = [0.86 0.86 0.86]; % subject lines 

% old RD 
% peak_T1 = squeeze(peakDataAve(:,1,:)); 
% peak_T2 = squeeze(peakDataAve(:,2,:));

% new KT 
sessionIdx = find(strcmp(pks.pkHeaders,'session')); 
peakDataAve = []; 
peakDataAve = mean(pks.pkVals,sessionIdx);
peak_T1 = squeeze(peakDataAve(:,:,:,1))'; 
peak_T2 = squeeze(peakDataAve(:,:,:,2))';

peakDiffSte(1) = std(peak_T1(1,idxPeaks)-peak_T1(2,idxPeaks),'omitnan')/sqrt(numel(idxPeaks));
peakDiffSte(2) = std(peak_T2(1,idxPeaks)-peak_T2(2,idxPeaks),'omitnan')/sqrt(numel(idxPeaks));

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
ylabel(sprintf('Normalized %s',variable))

% lines between T1 and T2 
if plotSubjects
    % subject lines 
    for i = idxPeaks
        s = plot([repmat(1+xJitter,[9,1]) repmat(2-xJitter,[9,1])],[peak_T1(2,i) peak_T2(1,i)],'Color',subjectColor);
    end
end
