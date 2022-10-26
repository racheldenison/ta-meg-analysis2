% kt_manuscriptFigs8_singleTrialPowerRevised.m

% September 28, 2022
% Single trial power supplemental
% Permutation test 

% Load singleTrialPower_all_s20_n10.mat

%% single trial amp per channel > ...
% square > nanmean channels > nanmean trials > ...
% nanmean sessions > nanmean subjects

% --- Data ---
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise');
% includeIdx = [1,2,3,4,6,7,8,9,10]; 
includeIdx = 1:10; 
saveFigs = 1; 

% condNames = {'cueT1','cueT2'};
condNames = {'all'};

normalize = 1; 

%% 
% --- Figure ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 
t = p.t(1):p.t(end); 
tUnit = 10; 
femtoUnit = 10e-15; 

clear groupG3 
% 1. single trial amplitudes per channel, Square 
for i = 1:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        groupG3.(cond).tfPowsTrials(:,i,:,:) = groupA(i).(cond).tfPowsTrials;
    end
end
% 2. Nanmean channels 
for iC = 1:numel(condNames)
    cond = condNames{iC};
    groupG3.(cond).tfPowsMeanCh = squeeze(mean(groupG3.(cond).tfPowsTrials,4,'omitnan'));
end
% 3. Nanmean trials 
for iC = 1:numel(condNames)
    cond = condNames{iC};
    groupG3.(cond).tfPowsMeanTrial = squeeze(mean(groupG3.(cond).tfPowsMeanCh,1,'omitnan'));
end
% 4. Nanmean sessions 
count = 1;
for i = 1:2:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        s1 = groupG3.(cond).tfPowsMeanTrial(i,:);
        s2 = groupG3.(cond).tfPowsMeanTrial(i+1,:);
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).powSubject(:,count) = s12mean;
    end
    count = count+1;
end
% 5. Select subjects
for iC = 1:numel(condNames)
    cond = condNames{iC};
    if normalize
        groupG3.(cond).powSubjectsSelected = groupG3.(cond).normSubject(:,includeIdx);
    else
        groupG3.(cond).powSubjectsSelected = groupG3.(cond).powSubject(:,includeIdx);
    end
end
% 6. Nanmean subjects 
for iC = 1:numel(condNames)
    cond = condNames{iC};
    groupG3.(cond).powGroup = mean(groupG3.(cond).powSubjectsSelected,2,'omitnan'); 
end

% Units 
val = groupG3.(cond).powGroup/femtoUnit^2; 

plot(t(1:tUnit:end),val(1:tUnit:end),'Color','k','LineWidth',2);

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle
% ytickformat('%.1e')
xlabel('Time (ms)')
ylabel('Power (fT^{2})')
xlim([-100 2400])

figTitle = sprintf('singleTrialPower_n%d_allTrials',numel(includeIdx)); 
if saveFigs
    saveas(gcf,sprintf('%s.svg', figTitle)) 
end

%% Check single subjects 
figure
set(gcf,'Position',[100 100 400 1200])
tUnit = 10;
for i = 1:numel(subjectNames)
    subplot (numel(subjectNames),1,i) 
    hold on 
    meg_figureStyle
    ylabel(und2space(subjectNames{i}))
    for iC = 1:numel(conds)
        cond = conds{iC};
        if normalize
            val = groupG3.(cond).normSubject(:,i);
        else
            val = groupG3.(cond).ITPCsubject(:,i);
        end
        l = plot(t(1:tUnit:end),val(1:tUnit:end),'Color',p.cueColors(iC,:),'LineWidth',2);
    end
end
xlabel('Time (ms)')

if normalize
    sgtitle('normalized ITPC')
    figTitle = sprintf('ITPC_byAtt_normalize_subjects');
else
    sgtitle('ITPC')
    figTitle = sprintf('ITPC_byAtt_subjects');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end

%% Do permutation test 
% Shuffle precue condition 1000x at subject level 

toi = [-100 2400];
tTOI = find(p.t==toi(1)):find(p.t==toi(end));

data1 = groupG3.cueT1.powSubjectsSelected(tTOI,:);
data2 = groupG3.cueT2.powSubjectsSelected(tTOI,:);

% clear unnecessary variables 
clear groupA, clear groupG, clear groupS

nPerm = 1000;
stat = 'sum'; 
variableName = 'singleTrialPower'; 

[A] = meg_nonparametricTest(data1,data2,nPerm,stat,variableName,tTOI);

%% 

