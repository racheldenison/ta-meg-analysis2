% 20 Hz ITPC 
% loading ITPC_20Hz_cue.mat

%% Settings 
condNames = 'precue'; % all
location = 'kantian'; % for file pathing, 'karen' for labmac, 'kantian' on laptop and iMac 

saveFigs = 0; 

normalize = 1; 
% baselineWin = 500:1000; % try old baseline window
baselineWin = 470:970;

% load R peak times from csv 
peakDat = sprintf('/Users/%s/Dropbox/Data/TANoise/fromRachel/itpcNorm_TS_Peaks_N10_20211225_workspace.mat',location); 
pk_RD = load(peakDat); 

% Read kt table peaks
peakDat = sprintf('/Users/%s/Dropbox/Data/TANoise/TANoise-stats/data/TANoise_peaks.csv',location);
pk_KT = readtable(peakDat); 

% Select peak data to use 
pk = pk_KT; 

% includeIdx = 1:10; % n = 10 (all) subjects to plot 
includeIdx = [1,2,3,4,6,7,8,9,10]; % n = 9 

%% Load data and directories 
switch condNames
    case 'precue'
        conds = {'cueT1','cueT2'}; 
        load('ITPC_20Hz_cue.mat') % for conds cue T1 cue T2 
    case 'all'
        conds = 'all';
end

p = meg_params('TANoise_ITPCsession8');
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise');
t = p.t(1):p.t(end); 
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/main/ITPC',location);

%% Compile group data 
% groupA = groupITPC_20Hz_avgTrial;
foi = 20; 
clear groupG3 
% 1. ITPC session
for i = 1:20
    for iC = 1:numel(conds)
        cond = conds{iC};
        groupG3.(cond).ITPCsession(:,i) = squeeze(A.(cond).session(foi,:,i));
    end
end
% 2. Normalize sessions 
if normalize
    baselineIdx = find(t==baselineWin(1)):find(t==baselineWin(end)); 
    for iC = 1:numel(conds)
    cond = conds{iC};
        val = groupG3.(cond).ITPCsession; 
        baselineVal = mean(val(baselineIdx,:),1,'omitnan');
        subBaseVal = val-baselineVal; 
        for iS = 1:numel(sessionNames)
            flipVal(:,iS) = subBaseVal(:,iS).*ITPCsession(iS);
        end
        groupG3.(cond).normSession = flipVal; 
    end
end
% 3. Nanmean sessions to subjects 
count = 1;
for iS = 1:2:20
    for iC = 1:numel(conds)
        cond = conds{iC};
        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).ITPCsession(:,iS)';
        s2 = groupG3.(cond).ITPCsession(:,iS+1)';
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).ITPCSubject(:,count) = s12mean;

        s1 = []; s2 = []; s12 = []; s12mean = []; 
        s1 = groupG3.(cond).normSession(:,iS)';
        s2 = groupG3.(cond).normSession(:,iS+1)';
        s12 = cat(1,s1,s2);
        s12mean = mean(s12,1,'omitnan');
        groupG3.(cond).normSubject(:,count) = s12mean;
    end
    count = count+1;
end
% 4. Select subjects
for iC = 1:numel(conds)
    cond = conds{iC};
    if normalize
        groupG3.(cond).ITPCSubjectsSelected = groupG3.(cond).normSubject(:,includeIdx);
    else
        groupG3.(cond).ITPCSubjectsSelected = groupG3.(cond).ITPCSubject(:,includeIdx);
    end  
end
% 5. Nanmean subjects 
for iC = 1:numel(conds)
    cond = conds{iC};
    groupG3.(cond).ITPCGroup = mean(groupG3.(cond).ITPCSubjectsSelected,2,'omitnan'); 
end

%% Figure (group) 
% --- Figure ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 

tUnit = 10; 

% Units 
val_cueT1 = groupG3.cueT1.ITPCGroup; 
val_cueT2 = groupG3.cueT2.ITPCGroup; 

plot(t(1:tUnit:end),val_cueT1(1:tUnit:end),'Color',p.cueColors(1,:),'LineWidth',2);
plot(t(1:tUnit:end),val_cueT2(1:tUnit:end),'Color',p.cueColors(2,:),'LineWidth',2);

% --- Format --- 
for i = 1:numel(p.eventTimes)
    xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
end
meg_figureStyle
% ytickformat('%.1e')
% yaxis to better match ITPC 
if normalize
    ylim([-0.04 0.1])
else
    ylim([0.25 0.4])
end
xlabel('Time (ms)')
ylabel('ITPC')
xlim([-100 2400])

if normalize
    figTitle = sprintf('ITPC_n%d_byAtt_normalize',numel(includeIdx));
else
    figTitle = sprintf('ITPC_n%d_byAtt',numel(includeIdx));
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

%% Peak ITPC times (subjects)
pkT = []; 
pkT(1,:) = [1207,... % 1 R0817
    1202,... % 2 R0898 
    1169,... % 3 R0959
    1205,... % 4 R0983
    NaN,... % 5 R1021
    1199,... % 6 R1103
    1175,... % 7 R1187
    1199,... % 8 R1373
    1116,... % 9 R1452
    1173]; % 10 R1507
pkT(2,:) = [1503,... % 1 R0817
    1512,... % 2 R0898 
    1477,... % 3 R0959
    1495,... % 4 R0983
    NaN,... % 5 R1021
    1449,... % 6 R1103
    1473,... % 7 R1187
    1502,... % 8 R1373
    1412,... % 9 R1452
    1484]; % 10 R1507

%% Figure (subjects) 
figure
set(gcf,'Position',[100 100 400 1200])
tUnit = 10;
sz = 80; 
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
            val = groupG3.(cond).ITPCSubject(:,i);
        end
        l = plot(t(1:tUnit:end),val(1:tUnit:end),'Color',p.cueColors(iC,:),'LineWidth',2);
        for iT = 1:2 % 2 targets
            pkTIdx = find(p.t==pkT(iT,i));
            scatter(pkT(iT,i),val(pkTIdx),sz,'^','MarkerFaceColor',p.cueColors(iC,:))         
        end
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

%% Setup permutation test 
% Shuffle precue condition 1000x at subject level 
normalize = 1;
dsample = 1; 
session = 0; % if not session, subject level s
includeIdx = [1 2 3 4 6 7 8 9 10]; % 1:10 

timeWin = 'fulltrial'; % 'targets', 'T1', 'T2', 'fulltrial'

switch timeWin
    case 'fulltrial'
        toi = [-100 2400];
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

data1 = []; data2 = []; 
if normalize
    variableName = 'ITPC_norm';  
    data1 = groupG3.cueT1.normSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.normSubject(tTOI,includeIdx);
    if session
        variableName = 'ITPC_norm_session';
        data1 = groupG3.cueT1.ITPCsession(tTOI,includeIdx); 
        data2 = groupG3.cueT2.ITPCsession(tTOI,includeIdx); 
    end
    if dsample
        tTOI = downsample(tTOI,10); 
        data1 = movmean(data1,10); 
        data2 = movmean(data2,10); 
        
        data1 = downsample(data1,10); 
        data2 = downsample(data2,10); 
    end
else
    variableName = 'ITPC'; 
    data1 = groupG3.cueT1.ITPCSubject(tTOI,includeIdx);
    data2 = groupG3.cueT2.ITPCSubject(tTOI,includeIdx);
end

%% Run permutation test 
clear groupA, clear groupG, clear groupS

nPerm = 1000;
stat = 'sum'; 
condName = 'precue'; 
ttestType = 'paired'; 
shuffleDim = 2;
figDim = 1; 
variable = 'att'; 

% [A,fH,figNames] = meg_permutationTest(data1,data2,nPerm,stat,ttestType, shuffleDim, figDim, variable);
[A] = meg_nonparametricTest(data1,data2,nPerm,stat,variableName,condName,tTOI);

%% Run peak analysis (sessions) 
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
                val = groupG3.(cond).ITPCSession(:,i+iS-1);
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
    sgtitle('normalized ITPC')
    figTitle = sprintf('ITPC_byAtt_norm_peaks_sessions');
else
    sgtitle('ITPC')
    figTitle = sprintf('ITPC_byAtt_peaks_sessions');
end
if saveFigs
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 
end

%% Sort peak data by RD subject idx
subjectNamesRD  = {'R0817',... % 1 
    'R1187',... % 2 
    'R0983',... % 3 
    'R0898',... % 4 
    'R1021',... % 5 
    'R1103',... % 6 
    'R0959',... % 7 
    'R1373',... % 8 
    'R1452',... % 9 
    'R1373'}; % 10
ktRDIdx = [1 4 7 3 5 6 2 8 9 10];
pkValsRD = NaN(10,2,2,2); 
for i = 1:numel(subjectNames)
    % save peaks RD idx 
    val = [];
    val = pkVals(ktRDIdx(i),:,:,:);
    disp(val)
    if ~isnan(val)
        pkValsRD(i,:,:,:) = val;
    end
end

%% For excel
iT = 2; iC = 2; iS = 2;
pkValsCSV = pkValsRD(:,iS,iC,iT);
disp([iT iC iS])


