% look at reaction time 
% April 15 2021 

%% load group behavioral data (n20)
load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/groupB.mat')

%% setup 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt);
p = meg_params('TANoise_ITPCsession8'); 

saveFigs = 1; 
if saveFigs 
    figDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/vis/figures/other'; 
    fileType = 'png'; % svg 
end

%% calculate reaction time mean and std 
for i = 1:numel(sessionNames)
    rt = groupB(i).rt; 
    rt_mean = nanmean(rt); 
    rt_std = nanstd(rt); 
    groupB(i).rt_mean = rt_mean; 
    groupB(i).rt_std = rt_std; 
end

%% average to subjects 
count = 0; 
for i = 1:2:numel(sessionNames)
    session1mean = groupB(i).rt_mean; % session 1 rt mean 
    session2mean = groupB(i+1).rt_mean; % session 2 rt mean 
    subjectMeanRt = mean([session1mean session2mean]); % subject rt mean 
    
    session1std = groupB(i).rt_std; % session 1 rt std
    session2std = groupB(i+1).rt_std; % session 2 rt std
    subjectStdRt = mean([session1std session2std]); % subject rt std 
    
    count = count + 1; 
    rtMeanSubjects(count) = subjectMeanRt; 
    rtStdSubjects(count) = subjectStdRt; 
end

%% separate rt by uppers and downers 
rtUpMean = []; 
rtUpStd = []; 
rtDownMean = []; 
rtDownStd = []; 
for i = 1:numel(subjectNames)
    if ITPCsubject(i)==1
        rtUpMean = [rtUpMean rtMeanSubjects(i)]; 
        rtUpStd = [rtUpStd rtStdSubjects(i)]; 
    elseif ITPCsubject(i)==-1
        rtDownMean = [rtDownMean rtMeanSubjects(i)]; 
        rtDownStd = [rtDownStd rtStdSubjects(i)]; 
    end
end  

%% choose variables to plot
variable = 'logrt'; % 'logrt'
disp(variable)
switch variable 
    case 'rt'
        % uppers 
        meanVar_up = rtUpMean; 
        stdVar_up = rtUpStd; 
        % downers
        meanVar_down = rtDownMean; 
        stdVar_down = rtDownStd; 
    case 'logrt'
        % uppers 
        meanVar_down = log(rtUpMean); 
        stdVar_down = log(rtUpStd); 
        % downers
        meanVar_down = log(rtDownMean); 
        stdVar_down = log(rtDownStd); 
    otherwise 
        error('variable not recognized (rt, logrt)');
end

%% plot 
figure % swarm charts of single trial rt's colored by uppers and downers 
set(gcf, 'Position',  [100, 100, 500, 300])

hold on 
for i = 1:numel(sessionNames)
    if ITPCsession(i)==1
        color = p.validColor;
    elseif ITPCsession(i)==-1
        color = p.invalidColor;
    elseif ITPCsession(i)==0
        color = [0.7 0.7 0.7]; 
    end
    rtFilter = groupB(i).rt(groupB(i).rt<5 & groupB(i).rt>0); 
    switch variable 
        case 'rt'
            swarmchart(repmat(i,size(rtFilter)),rtFilter,5,color,'MarkerEdgeAlpha',0.5,'MarkerFaceColor',color,'MarkerFaceAlpha',0.1);
        case 'logrt'
            swarmchart(repmat(i,size(rtFilter)),log(rtFilter),5,color,'MarkerEdgeAlpha',0.5,'MarkerFaceColor',color,'MarkerFaceAlpha',0.1);
            ylim([-2.5 1])
    end
end
ylabel('Reaction time (s)')
xlabel('Session')
xlim([0 numel(sessionNames)+1])
meg_figureStyle

if saveFigs 
    saveas(gcf,sprintf('%s/rt_sessions_distributions.%s',figDir,fileType)) 
end

%% plot mean and ste rt by uppers and dowwners 

markerSize = 14; 
lineWidth = 2; 
figure % rt mean and ste by uppers and downers 
set(gcf, 'Position',  [100, 100, 500, 300])
subplot 121 
hold on 
[h_rtmean,p_rtmean,ci_rtmean,stats_rtmean] = ttest2(meanVar_up,meanVar_down); % not significant 
[h_rtstd,p_rtstd,ci_rtstd,stats_rtstd] = ttest2(stdVar_up,stdVar_down); % not significant 
meg_figureStyle
plot(1, meanVar_up,  'o','MarkerSize',markerSize,'MarkerFaceColor',p.validColor,'MarkerEdgeColor','w','LineWidth',lineWidth) 
plot(2, meanVar_down,'o','MarkerSize',markerSize,'MarkerFaceColor',p.invalidColor,'MarkerEdgeColor','w','LineWidth',lineWidth)
xlim([0.5 2.5]) 
xticks([1 2])
xticklabels({'Uppers','Downers'})
ylabel('Mean reaction time (s)')
e = errorbar(1,mean(meanVar_up,'omitnan'),std(meanVar_up/sqrt(length(meanVar_up)),'omitnan'));
e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'k';
e.CapSize = 40;
e.LineWidth = 2; 
e = errorbar(2,mean(meanVar_down,'omitnan'),std(meanVar_down/sqrt(length(meanVar_down)),'omitnan'));
e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'k';
e.CapSize = 40;
e.LineWidth = 2; 

subplot 122 % rt std 
hold on 
meg_figureStyle
plot(1, stdVar_up,'o','MarkerSize',markerSize,'MarkerFaceColor',p.validColor,'MarkerEdgeColor','w','LineWidth',lineWidth) 
plot(2, stdVar_down,'o','MarkerSize',markerSize,'MarkerFaceColor',p.invalidColor,'MarkerEdgeColor','w','LineWidth',lineWidth)
xlim([0.5 2.5]) 
xticks([1 2])
xticklabels({'Uppers','Downers'})
ylabel('Std reaction time (s)')
e = errorbar(1,mean(stdVar_up,'omitnan'),std(stdVar_up/sqrt(length(stdVar_up)),'omitnan'));
e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'k';
e.CapSize = 40;
e.LineWidth = 2; 
e = errorbar(2,mean(stdVar_down,'omitnan'),std(stdVar_down/sqrt(length(stdVar_down)),'omitnan'));
e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'k';
e.CapSize = 40;
e.LineWidth = 2; 

if saveFigs 
    saveas(gcf,sprintf('%s/rt_mean_ste.%s',figDir,fileType)) 
end
