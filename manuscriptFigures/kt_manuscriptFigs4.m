% TANoise manuscript figs
% scatterplots 
% Behav and slope

%% Setup 
% Read data 
dataFile = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data/forMatlabFigs.xlsx'; 
T = readtable(dataFile); % subject level data 

% Figure settings 
xJitter = 0.2;
saveFigs = 0; 
subjectColor = [0.86 0.86 0.86]; 
figDir = '/Users/kantian/Dropbox/Data/TANoise/Manuscript_Figs/221004_updates'; 
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise'); 
capSize = 12; % errorbar cap size 
sSize = 30; % 20 subject scatter dot size 
gSize = 40; % group average scatter dot size 

%% Set variables
variableNames = string(T.Properties.VariableNames); 
idx = find(variableNames == "subject");
subject = T(:,idx); 

%% Behav: dprime 
dprime.T1.cueT1 = table2array(T(:,variableNames == "dprime_T1_C1")); 
dprime.T1.cueT2 = table2array(T(:,variableNames == "dprime_T1_C2")); 
dprime.T2.cueT1 = table2array(T(:,variableNames == "dprime_T2_C1")); 
dprime.T2.cueT2 = table2array(T(:,variableNames == "dprime_T2_C2")); 

% dprime figure 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
% Subject lines
for i = 1:numel(subjectNames)
    s = plot([1-xJitter 1+xJitter],[dprime.T1.cueT1 dprime.T1.cueT2],'Color',subjectColor); 
end
% T1 
% valid 
x = repmat(1-xJitter,[10,1]); 
e = errorbar(x(1),mean(dprime.T1.cueT1),std(dprime.T1.cueT1)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(7,:),'MarkerEdgeColor',p.cueColors(7,:),...
    'Color',p.cueColors(7,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,dprime.T1.cueT1,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,dprime.T1.cueT1,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid 
x = repmat(1+xJitter,[10,1]); 
e = errorbar(x(1),mean(dprime.T1.cueT2),std(dprime.T1.cueT2)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(8,:),'MarkerEdgeColor',p.cueColors(8,:),...
    'Color',p.cueColors(8,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,dprime.T1.cueT2,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,dprime.T1.cueT2,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% T2
% Subject lines
for i = 1:numel(subjectNames)
    s = plot([2-xJitter 2+xJitter],[dprime.T2.cueT2 dprime.T2.cueT1],'Color',subjectColor); 
end
% valid 
x = repmat(2-xJitter,[10,1]); 
e = errorbar(x(1),mean(dprime.T2.cueT2),std(dprime.T2.cueT2)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(7,:),'MarkerEdgeColor',p.cueColors(7,:),...
    'Color',p.cueColors(7,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,dprime.T2.cueT2,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,dprime.T2.cueT2,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid 
x = repmat(2+xJitter,[10,1]); 
e = errorbar(x(1),mean(dprime.T2.cueT1),std(dprime.T2.cueT1)/sqrt(10),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(8,:),'MarkerEdgeColor',p.cueColors(8,:),...
    'Color',p.cueColors(8,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,dprime.T2.cueT1,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,dprime.T2.cueT1,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
xlim([0 3])
ylabel('Sensitivity (d'')') 
xticks([1 2])
xticklabels({'T1','T2'})

if saveFigs
    figTitle = 'dPrime';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Behav: rt
rt.T1.cueT1 = table2array(T(:,variableNames == "rt_T1_C1"));
rt.T1.cueT2 = table2array(T(:,variableNames == "rt_T1_C2")); 
rt.T2.cueT1 = table2array(T(:,variableNames == "rt_T2_C1")); 
rt.T2.cueT2 = table2array(T(:,variableNames == "rt_T2_C2")); 

% dprime figure 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
% T1 
% Subject lines
for i = 1:numel(subjectNames)
    s = plot([1-xJitter 1+xJitter],[rt.T1.cueT1 rt.T1.cueT2],'Color',subjectColor); 
end
% valid 
x = repmat(1-xJitter,[10,1]); 
e = errorbar(x(1),mean(rt.T1.cueT1),std(rt.T1.cueT1)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(7,:),'MarkerEdgeColor',p.cueColors(7,:),...
    'Color',p.cueColors(7,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,rt.T1.cueT1,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,rt.T1.cueT1,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid 
x = repmat(1+xJitter,[10,1]); 
e = errorbar(x(1),mean(rt.T1.cueT2),std(rt.T1.cueT2)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(8,:),'MarkerEdgeColor',p.cueColors(8,:),...
    'Color',p.cueColors(8,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,rt.T1.cueT2,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,rt.T1.cueT2,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% T2
for i = 1:numel(subjectNames)
    s = plot([2-xJitter 2+xJitter],[rt.T2.cueT2 rt.T2.cueT1],'Color',subjectColor); 
end
% valid 
x = repmat(2-xJitter,[10,1]); 
e = errorbar(x(1),mean(rt.T2.cueT2),std(rt.T2.cueT2)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(7,:),'MarkerEdgeColor',p.cueColors(7,:),...
    'Color',p.cueColors(7,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,rt.T2.cueT2,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,rt.T2.cueT2,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid 
x = repmat(2+xJitter,[10,1]); 
e = errorbar(x(1),mean(rt.T2.cueT1),std(rt.T2.cueT1)/sqrt(10),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',p.cueColors(8,:),'MarkerEdgeColor',p.cueColors(8,:),...
    'Color',p.cueColors(8,:),'LineWidth',2);
e.CapSize = capSize; 
scatter(x,rt.T2.cueT1,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,rt.T2.cueT1,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

xlim([0 3])
ylabel('Reaction time (s)') 
xticks([1 2])
xticklabels({'T1','T2'})
ylim([0 1.5])

if saveFigs
    figTitle = 'rt';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Temporal expectation: slope, all trials 
var = table2array(T(:,variableNames == "slopes_allTrials"))*1000;
% var = slopes_allTrials*1000; % slopes (delta ITPC/s) 
xJitter = 0.12-0.05:0.05:0.12+0.05; 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
x = ones([numel(var),1]); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'Color','k','LineWidth',2);
e.CapSize = capSize; 
for i = 1:numel(var)
    xJit = randperm(numel(xJitter)); 
    scatter(x+xJitter(xJit(1)),var(i),sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',1,'MarkerEdgeColor','w','LineWidth',1)
end
yline(0,':k');
xticks([1])
xlim([0 2])
xticklabels({'All trials'})
ylabel('Slopes (delta ITPC/s)')
ylim([-0.02 0.14])
if saveFigs
    figTitle = 'slopes_allTrials';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Slope, by attention 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle

% Subject lines 
for i = 1:numel(subjectNames)
    var1 = table2array(T(:,variableNames == "slopes_cueT1"))*1000;
    var2 = table2array(T(:,variableNames == "slopes_cueT2"))*1000;
    s = plot([0 3],[var1 var2],'Color',subjectColor); 
end

% Precue T1 
var = table2array(T(:,variableNames == "slopes_cueT1"))*1000; % slopes (delta ITPC/s) 
x = repmat(0,[numel(var),1]); 
color = p.cueColors(1,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% Precue T2
var = table2array(T(:,variableNames == "slopes_cueT2"))*1000; % slopes (delta ITPC/s) 
x = repmat(3,[numel(var),1]); 
color = p.cueColors(2,:); 
errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

xlim([-1 4])
xticks([0 3])
xticklabels({'Precue T1', 'Precue T2'})
ylabel('Slopes (delta ITPC/s)')
% ylim([-0.02 0.14])
yline(0,':k')
if saveFigs
    figTitle = 'slopes_byCue';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Slope all trials, upper downer 
% Jan 11: need to get upper downer data!! insert 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
sSize = 40; 

% Upper 
var = table2array(T(:,variableNames == "slopes_allTrials_uppers"))*1000; % slopes (delta ITPC/s) 
var = var(~isnan(var)); 
x = repmat(0,[numel(var),1]); 
color = 'k'; 
errorbar(x(1),mean(var,'omitnan'),std(var,'omitnan')/sqrt(numel(var)),'Marker','o','MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
xJitter = -0.2:0.1:0.1; 
for i = 1:numel(var)
    xJit = randperm(numel(xJitter)); 
    scatter(x+xJitter(xJit(1)),var(i),sSize,'filled','MarkerFaceColor','w','MarkerEdgeColor','k', 'MarkerFaceAlpha',1)
end

% Downer 
var = table2array(T(:,variableNames == "slopes_allTrials_downers"))*1000; % slopes (delta ITPC/s) 
var = var(~isnan(var)); 
x = repmat(3,[numel(var),1]); 
color = 'k'; 
errorbar(x(1),mean(var,'omitnan'),std(var,'omitnan')/sqrt(numel(var)),'Marker','.','MarkerSize',40,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
scatter(x,var,sSize,'filled','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',1)
xlim([-1 4])
xticks([0 3])
xticklabels({'Upward', 'Downward'})
ylabel('Slopes (delta ITPC/s)')
% ylim([-0.02 0.14])
yline(0,':k')
if saveFigs
    figTitle = 'slopes_allTrials_upDown';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Slope, by jitter 
figure
hold on 
set(gcf,'Position',[100 100 400 400])
meg_figureStyle

% Subject lines 
for i = 1:numel(subjectNames)
    var1 = table2array(T(:,variableNames == "slopes_500"))*1000;
    var2 = table2array(T(:,variableNames == "slopes_700"))*1000;
    var3 = table2array(T(:,variableNames == "slopes_900"))*1000;
    var4 = table2array(T(:,variableNames == "slopes_1100"))*1000;
    var5 = table2array(T(:,variableNames == "slopes_1300"))*1000;
    var6 = table2array(T(:,variableNames == "slopes_1500"))*1000;
    s = plot([1:6],[var1 var2 var3 var4 var5 var6],'Color',subjectColor); 
end

% jitter 500
var = table2array(T(:,variableNames == "slopes_500"))*1000; % slopes (delta ITPC/s) 
x = repmat(1,[numel(var),1]); 
color = p.cueColors(1,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% 700
var = table2array(T(:,variableNames == "slopes_700"))*1000; % slopes (delta ITPC/s) 
x = repmat(2,[numel(var),1]); 
color = p.cueColors(2,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% 900
var = table2array(T(:,variableNames == "slopes_900"))*1000; % slopes (delta ITPC/s) 
x = repmat(3,[numel(var),1]); 
color = p.cueColors(3,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% 1100
var = table2array(T(:,variableNames == "slopes_1100"))*1000; % slopes (delta ITPC/s) 
x = repmat(4,[numel(var),1]); 
color = p.cueColors(4,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% 1300
var = table2array(T(:,variableNames == "slopes_1300"))*1000; % slopes (delta ITPC/s) 
x = repmat(5,[numel(var),1]); 
color = p.cueColors(5,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% 1500
var = table2array(T(:,variableNames == "slopes_1500"))*1000; % slopes (delta ITPC/s) 
x = repmat(6,[numel(var),1]); 
color = p.cueColors(6,:); 
e = errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,var,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

xlim([0 7])
xticks(1:6)
xticklabels({'500','700','900','1100','1300','1500'})
ylabel('Slopes (delta ITPC/s)')
% ylim([-0.02 0.14])
yline(0,':k')
if saveFigs
    figTitle = 'slopes_byJitter';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Behav: dprime by upper downer 
dprime.T1.cueT1 = table2array(T(:,variableNames == "dprime_T1_C1")); 
dprime.T1.cueT2 = table2array(T(:,variableNames == "dprime_T1_C2")); 
dprime.T2.cueT1 = table2array(T(:,variableNames == "dprime_T2_C1")); 
dprime.T2.cueT2 = table2array(T(:,variableNames == "dprime_T2_C2")); 

% dprime figure 
figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
gSize = 8; 
% Subject lines
% for i = 1:numel(subjectNames)
%     s = plot([1-xJitter 2-xJitter],[dprime.T1.cueT1(ITPCsubject==1) dprime.T1.cueT2(ITPCsubject==1)],'Color',subjectColor); 
% end
% T1 ----------
% valid (upper) 
color = p.cueColors(7,:); 
val = []; val = dprime.T1.cueT1(ITPCsubject==1); 
x = repmat(1-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% valid (downer) 
val = []; val = dprime.T1.cueT1(ITPCsubject==-1); 
x = repmat(1+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% invalid (upper) 
color = p.cueColors(8,:); 
val = []; val = dprime.T1.cueT2(ITPCsubject==1); 
x = repmat(2-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid (downer) 
val = []; val = dprime.T1.cueT2(ITPCsubject==-1); 
x = repmat(2+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% T2 ----------
% valid (upper) 
color = p.cueColors(7,:); 
val = []; val = dprime.T2.cueT2(ITPCsubject==1); 
x = repmat(4-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% valid (downer) 
val = []; val = dprime.T2.cueT2(ITPCsubject==-1); 
x = repmat(4+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% invalid (upper) 
color = p.cueColors(8,:); 
val = []; val = dprime.T2.cueT1(ITPCsubject==1); 
x = repmat(5-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid (downer) 
val = []; val = dprime.T2.cueT1(ITPCsubject==-1); 
x = repmat(5+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
xlim([0 6])
ylabel('Sensitivity (d'')') 
xticks([1.5 4.5])
xticklabels({'T1','T2'})

if saveFigs
    figTitle = 'dPrime_upDown';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Behav RT by up down

rt.T1.cueT1 = table2array(T(:,variableNames == "rt_T1_C1")); 
rt.T1.cueT2 = table2array(T(:,variableNames == "rt_T1_C2")); 
rt.T2.cueT1 = table2array(T(:,variableNames == "rt_T2_C1")); 
rt.T2.cueT2 = table2array(T(:,variableNames == "rt_T2_C2")); 

figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle
gSize = 8; 
% Subject lines
% for i = 1:numel(subjectNames)
%     s = plot([1-xJitter 2-xJitter],[rt.T1.cueT1(ITPCsubject==1) rt.T1.cueT2(ITPCsubject==1)],'Color',subjectColor); 
% end
% T1 ----------
% valid (upper) 
color = p.cueColors(7,:); 
val = []; val = rt.T1.cueT1(ITPCsubject==1); 
x = repmat(1-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% valid (downer) 
val = []; val = rt.T1.cueT1(ITPCsubject==-1); 
x = repmat(1+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% invalid (upper) 
color = p.cueColors(8,:); 
val = []; val = rt.T1.cueT2(ITPCsubject==1); 
x = repmat(2-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid (downer) 
val = []; val = rt.T1.cueT2(ITPCsubject==-1); 
x = repmat(2+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% T2 ----------
% valid (upper) 
color = p.cueColors(7,:); 
val = []; val = rt.T2.cueT2(ITPCsubject==1); 
x = repmat(4-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% valid (downer) 
val = []; val = rt.T2.cueT2(ITPCsubject==-1); 
x = repmat(4+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

% invalid (upper) 
color = p.cueColors(8,:); 
val = []; val = rt.T2.cueT1(ITPCsubject==1); 
x = repmat(5-xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,sSize,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')
% invalid (downer) 
val = []; val = rt.T2.cueT1(ITPCsubject==-1); 
x = repmat(5+xJitter,size(val)); 
e = errorbar(x(1),mean(val),std(val)/sqrt(size(val,1)),'Marker','o','MarkerSize',gSize,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
e.CapSize = capSize; 
scatter(x,val,'filled','MarkerFaceColor','w') % plot white first
scatter(x,val,sSize,'filled','MarkerFaceColor',color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w')

xlim([0 6])
yticks(0:0.5:1.5)
ylim([0 1.6])
ylabel('Reaction time (s)') 
xticks([1.5 4.5])
xticklabels({'T1','T2'})

if saveFigs
    figTitle = 'rt_upDown';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% peak T2 by session 
% load peak data workspace 'fromRachel' folder 
% peakData and peakData_dims

vals = squeeze(mean(peakData,1,'omitnan'));  % average cue 
figure
set(gcf,'Position',[100 100 600 400])

% --- T1 ---
subplot 121
hold on 
meg_figureStyle
val = []; 
val = squeeze(vals(1,:,:)); 
x = ones([10 1]); 
% subject lines 
for i = 1:10
    plot([x(1)*1 x(1)*2],[val(i,1) val(i,2)],'Color',subjectColor);  
end
scatter(x*1,val(:,1),gSize,'o','filled') % session 1 
scatter(x*2,val(:,2),gSize,'o','filled') % session 2
xlim([0 3])
ylabel('ITPC peak magnitude') 
xticks([1 2])
xticklabels({'Session 1','Session 2'})
xtickangle(30)
title('T1')

% --- T2 ---
subplot 122
hold on 
meg_figureStyle
val = []; 
val = squeeze(vals(2,:,:)); 
x = ones([10 1]); 
% subject lines 
for i = 1:10
    plot([x(1)*4 x(1)*5],[val(i,1) val(i,2)],'Color',subjectColor);  
end
scatter(x*4,val(:,1),gSize,'o','filled') % session 1 
scatter(x*5,val(:,2),gSize,'o','filled') % session 2
xlim([3 6])
xticks([4 5])
xticklabels({'Session 1','Session 2'})
xtickangle(30)
title('T2')
if saveFigs
    figTitle = 'peakSession';
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% baseline ITPC by upper downer 
% Read data 
dataFile = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data/TANoise_slopes_allTrials.csv'; 
T = readtable(dataFile); 
baseline = T.baselineITPC; 
count = 1; 
for i = 1:2:19
    val = [];
    val = mean(baseline(i:i+1)); 
    baselineS(count) = val; 
    count = count+1; 
end

figure
hold on 
set(gcf,'Position',[100 100 250 400])
meg_figureStyle

% Upper 
var = baselineS(ITPCsubject==1);  
x = repmat(0,[numel(var),1]); 
color = 'k'; 
errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','o','MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
xJitter = 0;
for i = 1:numel(var)
    xJit = randperm(numel(xJitter));
    scatter(x+xJitter(xJit(1)),var(i),sSize,'filled','MarkerFaceColor','w','MarkerEdgeColor','k', 'MarkerFaceAlpha',1)
end

% Downer 
var = baselineS(ITPCsubject==-1); 
x = repmat(3,[numel(var),1]); 
color = 'k'; 
errorbar(x(1),mean(var),std(var)/sqrt(numel(var)),'Marker','.','MarkerSize',40,'MarkerFaceColor',color,'MarkerEdgeColor',color,...
    'Color',color,'LineWidth',2);
scatter(x,var,sSize,'filled','MarkerFaceColor',color,'MarkerEdgeColor','w','MarkerFaceAlpha',1)

xlim([-1 4])
xticks([0 3])
xticklabels({'Upward', 'Downward'})
ylabel('Baseline ITPC')
% ylim([-0.02 0.14])
yline(0,':k')
if saveFigs
    figTitle = 'baselineITPC_upDown';
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Slopes old (Jan 19 2022) 
% count = 1; 
% for i = 1:2:19
%     val = [];
%     val = mean(T.slopesSessions(i:i+1)); 
%     slopesS(count) = val; 
%     
%     val = [];
%     val = mean(T.slopesSessionBuffer(i:i+1)); 
%     slopesBufferS(count) = val; 
%     
%     count = count+1; 
% end
% 
% x = ones(size(slopesS)); 
% figure
% hold on 
% scatter(x*1,slopesS)
% scatter(x*2,slopesBufferS) 
% meg_figureStyle
% xlim([0 3])
% 
% if saveFigs
%     figTitle = 'slopes_fit_oldnew';
%     saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
% end

%% Behavior (rt and dprime) by session (Jan 19 2022) 
% Read data 
dataFile = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data/TANoise_behav.csv'; 
T = readtable(dataFile); % session level data 

figure
subplot 221 
hold on 
sSize = 50;
session = 1; 
% ----- Session1 ----- 
%   --- T1 valid --- 
y = T.dPrime(T.Session==session & T.Target==1 & T.Cue==1); 
x = ones(size(y))*1; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 invalid --- 
y = T.dPrime(T.Session==session & T.Target==1 & T.Cue==2); 
x = ones(size(y))*2; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([1 2],[T.dPrime(T.Session==session & T.Target==1 & T.Cue==1) T.dPrime(T.Session==session & T.Target==1 & T.Cue==2)],'Color',subjectColor); 
end
%   --- T2 valid --- 
y = T.dPrime(T.Session==session & T.Target==2 & T.Cue==2); 
x = ones(size(y))*4; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 invalid --- 
y = T.dPrime(T.Session==session & T.Target==2 & T.Cue==1); 
x = ones(size(y))*5; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([4 5],[T.dPrime(T.Session==session & T.Target==2 & T.Cue==2) T.dPrime(T.Session==session & T.Target==2 & T.Cue==1)],'Color',subjectColor); 
end
%   --- styling --- 
meg_figureStyle
xticks([1.5,4.5])
xticklabels({'T1','T2'})
xlim([0 6])
ylim([0 4])
ylabel('Sensitivity (d'')')
title('Session 1') 

subplot 222
hold on 
session = 2; 
% ----- Session2 ----- 
%   --- T1 valid --- 
y = T.dPrime(T.Session==session & T.Target==1 & T.Cue==1); 
x = ones(size(y))*1; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 invalid --- 
y = T.dPrime(T.Session==session & T.Target==1 & T.Cue==2); 
x = ones(size(y))*2; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([1 2],[T.dPrime(T.Session==session & T.Target==1 & T.Cue==1) T.dPrime(T.Session==session & T.Target==1 & T.Cue==2)],'Color',subjectColor); 
end
%   --- T2 valid --- 
y = T.dPrime(T.Session==session & T.Target==2 & T.Cue==2); 
x = ones(size(y))*4; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 invalid --- 
y = T.dPrime(T.Session==session & T.Target==2 & T.Cue==1); 
x = ones(size(y))*5; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([4 5],[T.dPrime(T.Session==session & T.Target==2 & T.Cue==2) T.dPrime(T.Session==session & T.Target==2 & T.Cue==1)],'Color',subjectColor); 
end
%   --- styling --- 
meg_figureStyle
xticks([1.5,4.5])
xticklabels({'T1','T2'})
xlim([0 6])
ylim([0 4])
ylabel('Sensitivity (d'')')
title('Session 2')
lgd = legend({'valid','invalid'},'Box','on');

subplot 223
hold on 
sSize = 50;
session = 1; 
% ----- Session1 ----- 
%   --- T1 valid --- 
y = T.rt(T.Session==session & T.Target==1 & T.Cue==1); 
x = ones(size(y))*1; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 invalid --- 
y = T.rt(T.Session==session & T.Target==1 & T.Cue==2); 
x = ones(size(y))*2; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([1 2],[T.rt(T.Session==session & T.Target==1 & T.Cue==1) T.rt(T.Session==session & T.Target==1 & T.Cue==2)],'Color',subjectColor); 
end
%   --- T2 valid --- 
y = T.rt(T.Session==session & T.Target==2 & T.Cue==2); 
x = ones(size(y))*4; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 invalid --- 
y = T.rt(T.Session==session & T.Target==2 & T.Cue==1); 
x = ones(size(y))*5; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([4 5],[T.rt(T.Session==session & T.Target==2 & T.Cue==2) T.rt(T.Session==session & T.Target==2 & T.Cue==1)],'Color',subjectColor); 
end
%   --- styling --- 
meg_figureStyle
xticks([1.5,4.5])
xticklabels({'T1','T2'})
xlim([0 6])
ylabel('Reaction time (s)')
title('Session 1') 

subplot 224
hold on 
sSize = 50;
session = 2; 
% ----- Session1 ----- 
%   --- T1 valid --- 
y = T.rt(T.Session==session & T.Target==1 & T.Cue==1); 
x = ones(size(y))*1; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 invalid --- 
y = T.rt(T.Session==session & T.Target==1 & T.Cue==2); 
x = ones(size(y))*2; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T1 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([1 2],[T.rt(T.Session==session & T.Target==1 & T.Cue==1) T.rt(T.Session==session & T.Target==1 & T.Cue==2)],'Color',subjectColor); 
end
%   --- T2 valid --- 
y = T.rt(T.Session==session & T.Target==2 & T.Cue==2); 
x = ones(size(y))*4; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(7,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 invalid --- 
y = T.rt(T.Session==session & T.Target==2 & T.Cue==1); 
x = ones(size(y))*5; 
scatter(x,y,sSize,'filled','MarkerFaceColor',p.cueColors(8,:),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
%   --- T2 subject lines --- 
for i = 1:numel(subjectNames)
    s = plot([4 5],[T.rt(T.Session==session & T.Target==2 & T.Cue==2) T.rt(T.Session==session & T.Target==2 & T.Cue==1)],'Color',subjectColor); 
end
%   --- styling --- 
meg_figureStyle
xticks([1.5,4.5])
xticklabels({'T1','T2'})
xlim([0 6])
ylabel('Reaction time (s)')
title('Session 2') 

%% Slopes by fit time by session 
figure
subplot 211 
hold on 
meg_figureStyle
y_idx1 = 1:2:19; 
y_idx2 = 2:2:20; 

x = ones(size(y_idx1))*1; 
scatter(x*1,slopes(y_idx1)*1000,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
scatter(x*2,slopes(y_idx2)*1000,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 

scatter(x*4,slopes_buffer(y_idx1)*1000,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
scatter(x*5,slopes_buffer(y_idx2)*1000,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 

for i = 1:10
    s = plot([1 2 4 5],...
        [slopes(y_idx1(i))*1000 slopes(y_idx2(i))*1000 slopes_buffer(y_idx1(i))*1000 slopes_buffer(y_idx2(i))*1000],...
        'Color',subjectColor); 
end
xticks([1 2 4 5])
xticklabels({'Session 1','S2','S1','S2'})
xlim([0 6]) 
yline(0,':k') 
ylabel('Slope (delta ITPC/s)') 
title('No buffer                       Buffer') 

subplot 212 
hold on 
meg_figureStyle
slopes_noBuffer_vals(1,:) = slopes(y_idx1); 
slopes_noBuffer_vals(2,:) = slopes(y_idx2); 
y = mean(slopes_noBuffer_vals,1)*1000; 
x = ones(size(y)); 
scatter(x*1,y,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
slopes_buffer_vals(1,:) = slopes_buffer(y_idx1); 
slopes_buffer_vals(2,:) = slopes_buffer(y_idx2); 
y = mean(slopes_buffer_vals,1)*1000; 
scatter(x*2,y,sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
ylabel('Slope (delta ITPC/s)') 
xlim([0.5 2.5])
xticks([])
yline(0,':k') 
y1 = mean(slopes_noBuffer_vals,1)*1000;
y2 = mean(slopes_buffer_vals,1)*1000;
for i = 1:10
    s = plot([1 2],...
        [y1(i) y2(i)],...
        'Color',subjectColor);
end

if saveFigs
    figTitle = 'slopes_bufferCompare';
    saveas(gcf,sprintf('%s/interim/%s.png', figDir, figTitle))
end

%% 
dataFile = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data/TANoise_slopes_allTrials.csv'; 
T = readtable(dataFile); % subject level data 

val = []; 
count = 1; 
for i = 1:2:19 
    val1 = []; val2 = []; 
    val1 = mean(T.slopesSessionBufferS(i:i+1)); 
    val2 = mean(T.baselineITPC(i:i+1)); 
    val(1,count) = val1; % 1 precue to T1 slope 
    val(2,count) = val2; % baseline ITPC 
    count = count+1; 
end

figure
hold on 
x = ones(10,1); 
meg_figureStyle
% scatter(x*1,val(1,:),sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
% scatter(x*2,val(2,:),sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
scatter(val(1,:),val(2,:),sSize,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','w') 
xlabel('Precue to T1 slope (delta ITPC/s)') 
ylabel('Baseline ITPC') 
axis square
[rVal,pVal] = corrcoef(val(1,:),val(2,:));
xl = xlim;
yl = ylim; 
dim = [0.24 0.6 0.3 0.3];
str3 = 'n = 10'; 
str1 = sprintf('r = %0.3f',rVal(2));
str2 = sprintf('pVal = %0.3f',pVal(2));
str = {str3,str1,str2};
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','w');

if saveFigs
    figTitle = 'slopeVbaselineITPC';
    saveas(gcf,sprintf('%s/interim/%s.png', figDir, figTitle))
end
