function kt_manuscriptFigs7_singleTrialPower

%% Get data 
% run meg_runGroup
% change meg_runAnalysis settings with single trial ITPC 

% get groupA variable 
% get groupG variable (time x freq x subject) 
% singleTrialPower/singleTrialPower_n9.mat

%% Setup
condNames = {'cueT1','cueT2'};
condNames = {'all'};
saveFigs = 1; 

%% Average channels
for iS = 1:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        groupA(iS).(cond).meantfAmpsTrials = squeeze(mean(groupA(iS).(cond).tfPowsTrials,4,'omitnan')); 
        groupA(iS).(cond).meantfPowsTrials = squeeze(mean(groupA(iS).(cond).tfPowsTrials,4,'omitnan')); 
    end
end

%% First average trials test?? 
for iS = 1:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        groupA(iS).(cond).meantfAmpsTrials2 = squeeze(mean(groupA(iS).(cond).tfPowsTrials,1,'omitnan')); 
        groupA(iS).(cond).meantfPowsTrials2 = squeeze(mean(groupA(iS).(cond).tfPowsTrials,1,'omitnan')); 
    end
end

%% Average to subjects
clear groupS
count = 1;
averageFirst = 1;
averageTrialsFirst = 1; 
averageChannelsFirst = 0;
for iS = 1:2:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        if averageFirst
            if averageTrialsFirst
                % Average trials, then channels
                s1 = mean(groupA(iS).(cond).meantfPowsTrials2,2,'omitnan');
                s2 = mean(groupA(iS+1).(cond).meantfPowsTrials2,2,'omitnan');
                s12 = cat(2,s1,s2);
                groupS(count).(cond).tfPowsTrials = s12;
                groupS(count).(cond).meantfPowsTrials = squeeze(mean(s12,2,'omitnan'));
            elseif averageChannelsFirst
                % Average channels, then trials, per session
                s1 = mean(groupA(iS).(cond).meantfPowsTrials,1,'omitnan');
                s2 = mean(groupA(iS+1).(cond).meantfPowsTrials,1,'omitnan');
                s12 = cat(1,s1,s2);
                groupS(count).(cond).tfPowsTrials = s12;
                groupS(count).(cond).meantfPowsTrials = squeeze(mean(s12,1,'omitnan'));
            end
        else
            % Combine s1 s2 then average
            s1 = groupA(iS).(cond).meantfPowsTrials;
            s2 = groupA(iS+1).(cond).meantfPowsTrials;
            s12 = cat(1,s1,s2);
            groupS(count).(cond).tfPowsTrials = s12;
            groupS(count).(cond).meantfPowsTrials = squeeze(mean(s12,1,'omitnan'));
        end
    end
    count = count+1;
end

%% Make group matrix
clear groupG
for iS = 1:10
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        if iS==1
            groupG.(cond).meantfPowsTrials = [];
        end
        sVal = groupS(iS).(cond).meantfPowsTrials;
        groupG.(cond).meantfPowsTrials = cat(3,groupG.(cond).meantfPowsTrials,sVal);
    end
end

%% Plot 
% --- Data ---
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise');
% includeIdx = [1,2,3,4,6,7,8,9,10]; 
includeIdx = 1:10; 

% --- Figure ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 
t = p.t(1):p.t(end); 
tUnit = 10; 
femtoUnit = 10e-15; 
val_cueT1 = mean(groupG.cueT1.meantfPowsTrials(:,:,includeIdx),3,'omitnan'); 
val_cueT2 = mean(groupG.cueT2.meantfPowsTrials(:,:,includeIdx),3,'omitnan'); 
val_cueT1 = val_cueT1/femtoUnit^2; 
val_cueT2 = val_cueT2/femtoUnit^2; 
plot(t(1:tUnit:end),val_cueT1(1:tUnit:end),'Color',[p.cueColors(1,:)],'LineWidth',2);
plot(t(1:tUnit:end),val_cueT2(1:tUnit:end),'Color',[p.cueColors(2,:)],'LineWidth',2); 

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle
ytickformat('%.1e')
xlabel('Time (ms)')
ylabel('Power (fT^{2})')
% xlim([-100 2400])

figTitle = sprintf('singleTrialPower_n%d',numel(includeIdx)); 
if saveFigs
    saveas(gcf,sprintf('%s.svg', figTitle)) 
end

%% Plot ITPC by session 
clear groupG
s1Idx = 1:2:20; s1Idx(9)=[]; 
s2Idx = 2:2:20; s2Idx(10)=[]; 
condNames = {'all'}; 
for iS = 1:20
    for iC = 1:numel(condNames)
        cond = condNames{iC};
        if iS==1 
            groupG.(cond).ITPC = [];
        end
        sVal = groupA(iS).(cond).ITPCMean;
        groupG.(cond).ITPC = cat(2,groupG.(cond).ITPC,sVal);
    end
end

% Normalize
normIdx = find(p.t==471):find(p.t==970); 
baselines = groupG.all.ITPC(normIdx,:); 
for iS = 1:20
	baseline = mean(baselines(:,iS),'omitnan'); 
    normITPC = groupG.all.ITPC(:,iS)-baseline; 
    normITPCs = normITPC*ITPCsession(iS); 
    groupG.all.ITPCnorm(:,iS) = normITPCs; 
end

% --- Figure ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 
t = p.t(1):p.t(end); 
tUnit = 10; 

variable = 'ITPCnorm'; 

val_s1 = groupG.all.(variable)(:,s1Idx); 
val_s2 = groupG.all.(variable)(:,s2Idx); 


plot(t(1:tUnit:end),val_s1(1:tUnit:end,:),'Color',p.cueColors(6,:),'LineWidth',0.5);
plot(t(1:tUnit:end),val_s2(1:tUnit:end,:),'Color',p.cueColors(5,:),'LineWidth',0.5); 
% Means 
val_s1 = mean(val_s1,2,'omitnan');
val_s2 = mean(val_s2,2,'omitnan');
avg1 = plot(t(1:tUnit:end),val_s1(1:tUnit:end),'Color',p.cueColors(6,:),'LineWidth',2); % yellow 
avg2 = plot(t(1:tUnit:end),val_s2(1:tUnit:end),'Color',p.cueColors(5,:),'LineWidth',2); % green

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle
xlabel('Time (ms)')
ylabel('Normalized ITPC')
xlim([-100 2400])
legend([avg1 avg2],{'Session 1','Session 2'})

figTitle = sprintf('normalizedITPC_n%d',9); 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% ITPC peaks 
% filename = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data/TANoise_peaks.csv'; 
% T = readtable(filename);
% U = unstack(T,'PeakMag','Session');
% Anames = U.Properties.VariableNames([5 7 8 9]);
% A = table2array(U(:,[5 7 8 9]));
% 
% A = mean(A,2,'omitnan'); % collapse across cueing 
% 
% % Group fun
% varfun(@mean,T,'GroupingVariables',{'Cue','PeakMag'})

%% Single trial power plot subjects 
% load groupA by attention cond 

t = p.t(1):p.t(end); 
tUnit = 10; 
femtoUnit = 10e-15; 

for i = 1:numel(groupS)
    figure
    set(gcf,'Position',[100 100 600 400])
    hold on 
    subplot 311
    % Session
    val_T1 = groupA(sessionIdx(1)).cueT1.tfPows; 
    val_T1_ch = groupA(sessionIdx(1)).cueT1.tfPows; 
    plot(t(1:tUnit:end), val_T1(1:tUnit:end),'Color',[p.cueColors(1,:)],'LineWidth',2) 
    % All 
    val_T1 = groupA(sessionIdx(1)).cueT2.tfPows; 
    plot(t(1:tUnit:end), val_T1(1:tUnit:end),'Color',[p.cueColors(2,:)],'LineWidth',2) 
    xlim([-100 2400])
end

figTitle = sprintf('singleTrialPower_20Hz_att_%s',subjectNames{i}); 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% Compare single trial power R-K
% R0817_20171212 
% September 21, 2022
% load Rachel: analysis_singleTrials_R0817_TANoise_12.12.17_ebi_ft_topChannels5_allTrials_20Hz.mat
% load Karen: singleTrialPower_s20_n10.mat
xlims = [-1000,5000]; 

figure
sgtitle('Single trial power by precue')
subplot 211 % Rachel
hold on 
val = A.wAmpsAtt; 
val = squeeze(mean(val,2,'omitnan')); 
val1 = val(:,1); 
val2 = val(:,2); 
t = A.t;
plot(t,val1,'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,val2,'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 212 % Karen 
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPows;
val2 = groupA(sIdx).cueT2.tfPows; 
t = p.t; 
plot(t,val1,'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,val2,'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)
figTitle = 'R0817_20171212_SingleTrialPower_att'; 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% R-K check R0817_20171212 ITPC 
figure
sgtitle('ITPC by precue')
subplot 211 % Rachel
hold on 
val = A.wITPCAtt; 
val1 = val(:,1); 
val2 = val(:,2); 
t = A.t;
plot(t,val1,'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,val2,'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 212 % Karen 
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.ITPCMean;
val2 = groupA(sIdx).cueT2.ITPCMean; 
t = p.t; 
plot(t,val1,'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,val2,'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)
figTitle = 'R0817_20171212_ITPC_att'; 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% R-K check R0817_20171212
% single trial power, att, 5 channels,  average trials 
xlims = [-100 600];
figure
sgtitle('Single trial power by precue')
subplot 211 % Rachel
hold on 
val = A.wSpecAtt; % t x trials x att x ch
val = squeeze(mean(val,2,'omitnan')); 
val1 = val(:,1,:); 
val2 = val(:,2,:); 
t = A.t;
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',0.5)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',0.5)
meg_figureStyle
xlim(xlims)

subplot 212 % Karen 
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPowsTrials;
val2 = groupA(sIdx).cueT2.tfPowsTrials; 
val1 = mean(val1,1,'omitnan'); 
val2 = mean(val2,1,'omitnan'); 
t = p.t; 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',0.5)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',0.5)
meg_figureStyle
xlim(xlims)
figTitle = 'R0817_20171212_SingleTrialPower_att_ch'; 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% R-K check R0817_20171212
% single trial power, att, 1 channels, 3 trials 
xlims = [-1000,5000]; 
trialIdx = 1:5; 

figure
sgtitle(sprintf('Single trial power by precue\nTrials #%d-%d',trialIdx(1),trialIdx(end)))
subplot 211 % Rachel
hold on 
val = A.wAmpsAtt(:,trialIdx,:); % t x trials x att
val1 = val(:,:,1); 
val2 = val(:,:,2); 
t = A.t;
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',0.5)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',0.5)
meg_figureStyle
xlim(xlims)

subplot 212 % Karen 
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPowsTrials(trialIdx,:,:,:);
val2 = groupA(sIdx).cueT2.tfPowsTrials(trialIdx,:,:,:); 
val1 = mean(val1,4,'omitnan'); 
val2 = mean(val2,4,'omitnan'); 
t = p.t; 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',0.5)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',0.5)
meg_figureStyle
xlim(xlims)
figTitle = 'R0817_20171212_SingleTrialPower_att_ch'; 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% Compare order of averaging 
figure

subplot 311 % Karen, first average trials, then channels 
title('KT, first average trials, then channels')
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPowsTrials;
val2 = groupA(sIdx).cueT2.tfPowsTrials; 
val1 = mean(val1,1,'omitnan'); % average trials 
val2 = mean(val2,1,'omitnan'); 
val1 = mean(val1,4,'omitnan'); % average channels 
val2 = mean(val2,4,'omitnan'); 
t = p.t; 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 312 % Karen, first average channels, then trials 
title('KT, first average channels, then trials')
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPowsTrials;
val2 = groupA(sIdx).cueT2.tfPowsTrials; 
val1 = mean(val1,4,'omitnan'); % average channels 
val2 = mean(val2,4,'omitnan'); 
val1 = mean(val1,1,'omitnan'); % average trials 
val2 = mean(val2,1,'omitnan'); 
t = p.t; 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 313 % Karen, average all
title('KT, average all')
hold on 
sIdx = 1; % session index 
val1 = groupA(sIdx).cueT1.tfPowsTrials;
val2 = groupA(sIdx).cueT2.tfPowsTrials; 
val1 = mean(val1,'all','omitnan'); % average channels 
val2 = mean(val2,'all','omitnan'); 
t = p.t; 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',2)
meg_figureStyle
xlim(xlims)

figTitle = 'R0817_20171212_SingleTrialPower_att_compareAverageCh-averageTrialFirst'; 
if saveFigs
    saveas(gcf,sprintf('%s.png', figTitle)) 
end

%% Try other variables
figure

subplot 311
title('RD, wAmpsAll')
hold on 
t = A.t;
val = A.wAmpsAll; 
val = mean(val,2,'omitnan');
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 312
title('RD, wSpecAll')
hold on 
t = A.t;
val = A.wSpecAll; 
val = mean(val,2,'omitnan');
val = mean(val,3,'omitnan'); 
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlim(xlims)

subplot 313
title('RD, powAll = wAmpsAll^2')
hold on 
t = A.t;
val = A.wAmpsAll.^2; 
val = mean(val,2,'omitnan');
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlim(xlims)

%% 
figure
subplot 311
title('RD, stfAmpsAtt')
hold on 
t = A.t(1):10:A.t(end);
val = A.stfAmpsAtt; 
val1 = val(20,:,1); 
val2 = val(20,:,2); 
plot(t,squeeze(val1),'Color',p.cueColors(1,:),'LineWidth',2)
plot(t,squeeze(val2),'Color',p.cueColors(2,:),'LineWidth',2)

meg_figureStyle
xlim(xlims)


%% Save data matrices
dataDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/main/singleTrialPower';
dataName = 'singleTrialPower_s20_n10'; 
filename = sprintf('%s/%s.mat',dataDir,dataName); 
save(filename,'groupA','groupS','groupG','-v7.3')


