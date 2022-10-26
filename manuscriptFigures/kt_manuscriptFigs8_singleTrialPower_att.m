% kt_manuscriptFigs8_singleTrialPowerRevised.m

% September 26, 2022
% Single trial power supplemental
% Permutation test 

% Load singleTrialPower_att_s20_n10.mat

figDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/main/singleTrialPower';

%% single trial amp per channel > ...
% square > nanmean channels > nanmean trials > ...
% nanmean sessions > nanmean subjects

% --- Data ---
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise');
% includeIdx = [1,2,3,4,6,7,8,9,10]; 
includeIdx = 1:10; 
saveFigs = 1; 
normalize = 1; 

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
    groupG3.cueT1.tfPowsTrials(:,i,:,:) = groupA(i).cueT1.tfPowsTrials; 
    groupG3.cueT2.tfPowsTrials(:,i,:,:) = groupA(i).cueT2.tfPowsTrials; 
end
% 2. Nanmean channels 
groupG3.cueT1.tfPowsMeanCh = squeeze(mean(groupG3.cueT1.tfPowsTrials,4,'omitnan')); 
groupG3.cueT2.tfPowsMeanCh = squeeze(mean(groupG3.cueT2.tfPowsTrials,4,'omitnan')); 
% 3. Nanmean trials 
groupG3.cueT1.powSession = squeeze(mean(groupG3.cueT1.tfPowsMeanCh,1,'omitnan')); 
groupG3.cueT2.powSession = squeeze(mean(groupG3.cueT2.tfPowsMeanCh,1,'omitnan')); 
% Optional normalize 
if normalize
    baselineWin = 470:970; 
    baselineIdx = find(t==baselineWin(1)):find(t==baselineWin(end)); 
    for iC = 1:numel(conds)
    cond = conds{iC};
        val = groupG3.(cond).powSession; 
        baselineVal = mean(val(:,baselineIdx),2,'omitnan');
        subBaseVal = val-baselineVal; 
        for iS = 1:20
            flipVal(:,iS) = subBaseVal(iS,:).*ITPCsession(iS);
        end
        groupG3.(cond).normSession = flipVal'; 
    end
end
% 4. Nanmean sessions to subjects 
count = 1;
for iS = 1:2:20
    for iC = 1:numel(conds)
        cond = conds{iC};
        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).powSession(iS,:);
        s2 = groupG3.(cond).powSession(iS+1,:);
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).powSubject(:,count) = s12mean;

        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).normSession(iS,:);
        s2 = groupG3.(cond).normSession(iS+1,:);
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).normSubject(:,count) = s12mean;
    end
    count = count+1;
end
% 5. Select subjects
for iC = 1:numel(conds)
    cond = conds{iC};
    if normalize 
        groupG3.(cond).powSubjectsSelected = groupG3.(cond).normSubject(:,includeIdx); 
    else
        groupG3.(cond).powSubjectsSelected = groupG3.(cond).powSubject(:,includeIdx); 
    end
end
% 6. Nanmean subjects 
for iC = 1:numel(conds)
    cond = conds{iC};
    groupG3.(cond).powGroup = mean(groupG3.(cond).powSubjectsSelected,2,'omitnan'); 
end
% Units 
val_cueT1 = mean(groupG3.cueT1.normSubject(:,includeIdx),2,'omitnan')/femtoUnit^2; 
val_cueT2 = mean(groupG3.cueT2.normSubject(:,includeIdx),2,'omitnan')/femtoUnit^2; 

plot(t(1:tUnit:end),val_cueT1(1:tUnit:end),'Color',[p.cueColors(1,:)],'LineWidth',2);
plot(t(1:tUnit:end),val_cueT2(1:tUnit:end),'Color',[p.cueColors(2,:)],'LineWidth',2); 

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle
% ytickformat('%.1e')
xlabel('Time (ms)')
if normalize
    ylabel('Normalized power (fT^{2})')
else
    ylabel('Power (fT^{2})')
end
xlim([-100 2400])

figTitle = sprintf('singleTrialPower_n%d_normalized',numel(includeIdx)); 
if saveFigs
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

%% Single subject figures 
normalize = 1; 

figure
set(gcf,'Position',[100 100 400 1200])
tUnit = 10;
alpha = 0.7; 
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
            val = groupG3.(cond).powSubject(:,i);
        end
        l = plot(t(1:tUnit:end),val(1:tUnit:end)/femtoUnit^2,'Color',p.cueColors(iC,:),'LineWidth',2);
        l.Color(4) = alpha; 
    end
    xlim([t(1) t(end)])
end
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
xlabel('Time (ms)')

if normalize
    sgtitle('Normalized single trial power')
    figTitle = sprintf('singleTrialPower_byAtt_normalize_subjects');
else
    sgtitle('Single trial power')
    figTitle = sprintf('singleTrialPower_byAtt_subjects');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end

%% Do permutation test 
% Shuffle precue condition 1000x at subject level 

normalize = 1;

toi = [-100 2400];% redo restrict toi to T1:T2+300
% toi = p.eventTimes(2):p.eventTimes(3)+300; % restrict toi to T1:T2+300
tTOI = find(p.t==toi(1)):find(p.t==toi(end));

if normalize
    variableName = 'singleTrialPower_norm';  
    includeIdx = [1,2,3,4,6,7,8,9,10]; 
    data1 = groupG3.cueT1.normSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.normSubject(tTOI,includeIdx);
else
    variableName = 'singleTrialPower'; 
    includeIdx = 1:10;
    data1 = groupG3.cueT1.powSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.powSubject(tTOI,includeIdx);
end

% clear unnecessary variables 
clear groupA, clear groupG, clear groupS

nPerm = 1000;
stat = 'sum'; 
variableName = 'singleTrialPower'; 

condName = 'precue'; 

[A] = meg_nonparametricTest(data1,data2,nPerm,stat,variableName,condName,tTOI);

%% Run peak analysis (sessions) 
% load rd peak data structure 
clear pkVals 
load('/Users/kantian/Dropbox/github/ta-meg-analysis2/main/ITPC/peaks.mat') % for peak times ;
figure
set(gcf,'Position',[100 100 700 1200])
tUnit = 10;
sz = 80; 
pkVals = NaN(10,2,2,2); % subject x session (2) x cue (2) x target (2) 
for i = 1:2:numel(sessionNames) % Sessions (by subject counter)
    for iC = 1:numel(conds) % cue conds 
        cond = conds{iC};       
        for iS = 1:2 % sessions
            if normalize
                val = groupG3.(cond).normSession(i+iS-1,:);
            else
                val = groupG3.(cond).powSession(i+iS-1,:);
            end
            subplot (numel(subjectNames),2,i+iS-1)
            hold on
            meg_figureStyle
            xlabel(und2space(sessionNames{i}))
            
            l = plot(t(1:tUnit:end),val(1:tUnit:end),'Color',p.cueColors(iC,:),'LineWidth',2);
            for iT = 1:2 % 2 targets
                pkTIdx = find(p.t==pkT(iT,(i+1)/2));
                pkTIdx = pkTIdx-50:pkTIdx+50; % window around pk time to average 
                pkVal = mean(val(pkTIdx),'omitnan'); 
                scatter(pkT(iT,(i+1)/2),pkVal,sz,'^','MarkerFaceColor',p.cueColors(iC,:))
                % save peaks
                disp(pkVal) 
                if ~isnan(pkVal)
                    pkVals((i+1)/2,iS,iC,iT) = pkVal;
                end
            end
        end
        
    end
end
xlabel('Time (ms)')

pkHeaders = {'subject',...
    'session',...
    'cue',...
    'target'}; 
% save mat
save(sprintf('%s/peaks2.mat',figDir),'pkVals','pkT','pkHeaders')
if normalize
    sgtitle('normalized single trial power')
    figTitle = sprintf('singleTrialPow_byAtt_norm_peaks_sessions');
else
    sgtitle('ITPC')
    figTitle = sprintf('singleTrialPow_byAtt_peaks_sessions');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end





