
% load data groupITPC_20Hz_avgTrial (all trials) 
% load averageTrialPower_n10_cue.mat
% groupA = groupITPC_20Hz_avgTrial;

%% Settings
normalize = 1; 
baselineWin = 470:970; 

condName = 'att'; %'all' 'att'

saveFigs = 1; 

% includeIdx = 1:10; % subjects to plot 
includeIdx = [1,2,3,4,6,7,8,9,10]; 

%% Load data and directories 
p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise');
t = p.t(1):p.t(end); 

switch condName
    case 'att'
        conds = {'cueT1','cueT2'};
        load('averageTrialPower_s20_cue.mat')
    case 'all'
        conds = {'all'}; 
        load('averageTrialPower_s20_allTrials.mat') 
end
figDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/main/averageTrialPower';

%% Compile group data 
clear groupG3 
% 1. average trial amplitudes per channel, all trials, square 
for i = 1:20
    for iC = 1:numel(conds)
        cond = conds{iC};
        groupG3.(cond).tfPowsTrials(:,i,:) = groupA(i).(cond).tfPows;
    end
end
% 2. Nanmean channels 
for iC = 1:numel(conds)
    cond = conds{iC};
    groupG3.(cond).powSession = squeeze(mean(groupG3.(cond).tfPowsTrials,3,'omitnan'));
end
% 3. Optional normalize 
baselineIdx = find(t==baselineWin(1)):find(t==baselineWin(end));
for iC = 1:numel(conds)
    cond = conds{iC};
    val = groupG3.(cond).powSession;
    baselineVal = mean(val(baselineIdx,:),1,'omitnan');
    subBaseVal = val-baselineVal;
    for iS = 1:numel(sessionNames)
        flipVal(:,iS) = subBaseVal(:,iS).*ITPCsession(iS);
    end
    groupG3.(cond).normSession = flipVal;
end
% 4. Nanmean sessions to subjects 
count = 1;
for iS = 1:2:20
    for iC = 1:numel(conds)
        cond = conds{iC};
        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).powSession(:,iS)';
        s2 = groupG3.(cond).powSession(:,iS+1)';
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).powSubject(:,count) = s12mean;

        % normalized
        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).normSession(:,iS)';
        s2 = groupG3.(cond).normSession(:,iS+1)';
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).normSubject(:,count) = s12mean;
    end
    count = count+1;
end
% 5. Select subjects
for iC = 1:numel(conds)
    cond = conds{iC};
    groupG3.(cond).powSubjectsSelected = groupG3.(cond).powSubject(:,includeIdx); 
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

%% Figure (group) 
% --- Figure ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 

tUnit = 10; 
femtoUnit = 10e-15; 

switch condName
    case 'att'
        val_cueT1 = groupG3.cueT1.powGroup/femtoUnit^2; 
        val_cueT2 = groupG3.cueT2.powGroup/femtoUnit^2; 
        plot(t(1:tUnit:end),val_cueT1(1:tUnit:end),'Color',p.cueColors(1,:),'LineWidth',2);
        plot(t(1:tUnit:end),val_cueT2(1:tUnit:end),'Color',p.cueColors(2,:),'LineWidth',2);

        % scale axes to match ITPC 
        if normalize
            ylim([-300 700]) % yaxis to better match ITPC
        else
            ylim([500 1400])
            yticks([600:200:1400])
        end

        % significance bar (permutation corrected) 
        load('A_permutationTest_221007_n9_att_avgTrialPower.mat')
        tPermute = -100:2400; 
        sigTimesIdx = (find(A.h_true)); 
        sigTimes = tPermute(sigTimesIdx);
        yl = ylim; 
        plot(sigTimes,ones(size(sigTimes))*yl(2)*0.95,'k','LineWidth',2)

    case 'all'
        val_all = groupG3.all.powGroup/femtoUnit^2; 
        plot(t(1:tUnit:end),val_all(1:tUnit:end),'Color','k','LineWidth',2);
        ylim([500 1300])
end

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle

xlabel('Time (ms)')
ylabel('Power (fT^{2})')
xlim([-100 2400])

if normalize
    figTitle = sprintf('averageTrialPower_n%d_%s_normalize',numel(includeIdx),condName);
else
    figTitle = sprintf('averageTrialPower_n%d_%s',numel(includeIdx),condName);
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

%% Figure (subjects) 
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
            val = groupG3.(cond).powSubject(:,i);
        end
        if numel(conds)==1
            l = plot(t(1:tUnit:end),val(1:tUnit:end),'Color','k','LineWidth',2);
        else
            l = plot(t(1:tUnit:end),val(1:tUnit:end),'Color',p.cueColors(iC,:),'LineWidth',2);
        end
    end
end
xlabel('Time (ms)')

if normalize
    sgtitle('Normalized average trial power')
    figTitle = sprintf('avgTrialPower_byAtt_normalize_subjects');
else
    sgtitle('Average trial power')
    figTitle = sprintf('avgTrialPower_byAtt_subjects');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end

%% Setup permutation test 
% Shuffle precue condition 1000x at subject level 
normalize = 1;
session = 0; % if not session, subject level s
includeIdx = [1 2 3 4 6 7 8 9 10]; % 1:10 

timeWin = 'fulltrial'; % 'targets', 'T1', 'T2', 'fulltrial'

switch timeWin
    case 'fulltrial'
        toi = [p.eventTimes(1) p.eventTimes(4)]; % [0 2300]
    case 'targets'
        toi = p.eventTimes(2):p.eventTimes(3)+300; 
    case 'T1'
        toi = p.eventTimes(2):p.eventTimes(2)+300; 
    case 'T2'
        toi = p.eventTimes(3):p.eventTimes(3)+300; 
    otherwise
        error('Time window not specified') 
end
tTOI = find(p.t==toi(1)):find(p.t==toi(end));

if normalize
    variableName = 'avgTrialPower_norm';  
    data1 = groupG3.cueT1.normSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.normSubject(tTOI,includeIdx);
    if session
        variableName = 'avgTrialPower_norm_session';
        data1 = groupG3.cueT1.tfPows(tTOI,includeIdx); 
        data2 = groupG3.cueT2.tfPows(tTOI,includeIdx); 
    end
else
    variableName = 'avgTrialPower'; 
    data1 = groupG3.cueT1.powSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.powSubject(tTOI,includeIdx);
end

%% Do permutation test 

% clear unnecessary variables 
clear groupG, clear groupS

nPerm = 1000;
stat = 'sum'; 
variableName = 'averageTrialPower_norm'; 
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
                val = groupG3.(cond).normSession(:,i+iS-1);
            else
                val = groupG3.(cond).powSession(:,i+iS-1);
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
save(sprintf('%s/peaks_AverageTrial.mat',figDir),'pkVals','pkT','pkHeaders')
if normalize
    sgtitle('normalized average trial power')
    figTitle = sprintf('avgTrialPow_byAtt_norm_peaks_sessions');
else
    sgtitle('Average Trial Power')
    figTitle = sprintf('avgTrialPow_byAtt_peaks_sessions');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end
 


