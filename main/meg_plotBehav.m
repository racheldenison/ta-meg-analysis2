function meg_plotBehav(expt, user)

% function meg_plotBehav(expt, user)
% plots sensitivity, accuracy, and rt by condition by axis tilt (optional) 
%
% INPUTS
% expt
%   string, 'TA2' or 'TANoise'
% user (optional)
%   string giving the user name to determine path to save figures
%
% Karen Tian
% May 2020

%% inputs
expt = 'TA2'; % TANoise 
user = 'karen'; % for file path to save figures ß
plotAxis = 1; % slice by axis orientation 
saveFigs = 1; % saves .svg to current directory 

%% setup
[sessionNames,subjectNames] = meg_sessions(expt); 
targets = {'T1','T2'}; 
p = meg_params(sprintf('%s_Analysis',expt)); 

%% make groupDAll(session).D (data)...I (index)...B (behav) 
groupDAll = []; 
for i = 1:numel(sessionNames)
    [~, D, ~, I, B] = meg_runAnalysis(expt, sessionNames{i}, user); 
    disp(sessionNames{i})
    % make group 
    groupDAll(i).slice = 'cue'; 
    groupDAll(i).D = D; 
    groupDAll(i).I = I; 
    groupDAll(i).B = B; 
end

%% make groupD.field(time x ch x session) 

fields = fieldnames(D); 
groupD = []; 
for i = 1:numel(sessionNames)
    for iF = 1:numel(fields)
        groupD.(fields{iF})(:,:,i) = nanmean(groupDAll(i).D.(fields{iF}),3); 
    end
    disp(i)
end

%% make groupB.(target).(cue)

groupB = [];
for i = 1:numel(sessionNames)
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            targetIdx = find(groupDAll(i).B.responseTarget==iT);
            fieldIdx = groupDAll(i).I.(fields{iF}); 
            unionIdx = intersect(targetIdx,fieldIdx); 
            
            groupB.(targets{iT}).(fields{iF}).acc(:,i) = groupDAll(i).B.acc(unionIdx);
            groupB.(targets{iT}).(fields{iF}).rt(:,i) = groupDAll(i).B.rt(unionIdx);
            groupB.(targets{iT}).(fields{iF}).targetOrientation(:,i) = groupDAll(i).B.targetOrientation(unionIdx);
            
            groupB.(targets{iT}).(fields{iF}).hit(:,i) = groupDAll(i).B.discrimHMFC(unionIdx,1);
            groupB.(targets{iT}).(fields{iF}).fa(:,i)  = groupDAll(i).B.discrimHMFC(unionIdx,3);     
            
            if plotAxis
                % orientation
                groupB.(targets{iT}).(fields{iF}).axis(:,i) = groupDAll(i).B.t1t2Axes(unionIdx,iT);
                % accuracy
                acc = groupDAll(i).B.acc(unionIdx);
                groupB.(targets{iT}).(fields{iF}).acc0(:,i) = acc(groupB.(targets{iT}).(fields{iF}).axis(:,i)==0); % 0=vertical?
                groupB.(targets{iT}).(fields{iF}).acc90(:,i) = acc(groupB.(targets{iT}).(fields{iF}).axis(:,i)==90); % 90=horizontal?
                % rt
                rt = groupDAll(i).B.rt(unionIdx);
                groupB.(targets{iT}).(fields{iF}).rt0(:,i) = rt(groupB.(targets{iT}).(fields{iF}).axis(:,i)==0); % 0=vertical?
                groupB.(targets{iT}).(fields{iF}).rt90(:,i) = rt(groupB.(targets{iT}).(fields{iF}).axis(:,i)==90); % 90=horizontal?
                % hit
                hit = groupDAll(i).B.discrimHMFC(unionIdx,1);
                groupB.(targets{iT}).(fields{iF}).hit0(:,i) = hit(groupB.(targets{iT}).(fields{iF}).axis(:,i)==0); % 0=vertical?
                groupB.(targets{iT}).(fields{iF}).hit90(:,i) = hit(groupB.(targets{iT}).(fields{iF}).axis(:,i)==90); % 90=horizontal?
                % fa
                fa = groupDAll(i).B.discrimHMFC(unionIdx,3);
                groupB.(targets{iT}).(fields{iF}).fa0(:,i) = fa(groupB.(targets{iT}).(fields{iF}).axis(:,i)==0); % 0=vertical?
                groupB.(targets{iT}).(fields{iF}).fa90(:,i) = fa(groupB.(targets{iT}).(fields{iF}).axis(:,i)==90); % 90=horizontal?
            end
        end
    end
end

%% avg behavior sessions to subjects

for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        nTrialsYes = size(groupB.(targets{iT}).(fields{iF}).acc,1)/2; 
        countSubject = 1;
        for i = 1:2:numel(sessionNames)
            groupB.(targets{iT}).(fields{iF}).subjectAcc(:,countSubject) = nanmean(nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).acc(:,i:i+1),2)));
            groupB.(targets{iT}).(fields{iF}).subjectRT(:,countSubject) = nanmean(nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).rt(:,i:i+1),2)));
            
            % groupB.(targets{iT}).(fields{iF}).subjectHit(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).hit(:,i:i+1),2);
            % groupB.(targets{iT}).(fields{iF}).subjectFA(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).fa(:,i:i+1),2);
            
            % dprime 
            hit = nansum(groupB.(targets{iT}).(fields{iF}).hit)./nTrialsYes; 
            fa = nansum(groupB.(targets{iT}).(fields{iF}).fa)./nTrialsYes; 
            
            hit(hit==1) = 0.99; 
            hit(hit==0) = 0.01; 
            fa(fa==1) = 0.99; 
            fa(fa==0) = 0.01;
            
            groupB.(targets{iT}).(fields{iF}).sessionDprime = norminv(hit) - norminv(fa); 
            groupB.(targets{iT}).(fields{iF}).subjectDprime(:,countSubject) = nanmean(norminv(hit(:,i:i+1)) - norminv(fa(:,i:i+1)));
            groupB.(targets{iT}).(fields{iF}).dPrime = nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime);
            groupB.(targets{iT}).(fields{iF}).dPrimeSte = nanstd(groupB.(targets{iT}).(fields{iF}).subjectDprime)/sqrt(numel(subjectNames));
            
            if plotAxis 
                groupB.(targets{iT}).(fields{iF}).subjectAcc0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).acc0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectAcc90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).acc90(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectRT0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).rt0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectRT90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).rt90(:,i:i+1),2);
               
                % axis 0, vertical 
                hit0 = nansum(groupB.(targets{iT}).(fields{iF}).hit0)./(nTrialsYes/2); 
                fa0 = nansum(groupB.(targets{iT}).(fields{iF}).fa0)./(nTrialsYes/2); 
                
                hit0(hit0==1) = 0.99;
                hit0(hit0==0) = 0.01;
                fa0(fa0==1) = 0.99;
                fa0(fa0==0) = 0.01;
                groupB.(targets{iT}).(fields{iF}).sessionDprime0 = norminv(hit0)- norminv(fa0);
                groupB.(targets{iT}).(fields{iF}).subjectDprime0(:,countSubject) = nanmean(norminv(hit0(:,i:i+1)) - norminv(fa0(:,i:i+1)));
                
                % axis 90, horizontal 
                hit90 = nansum(groupB.(targets{iT}).(fields{iF}).hit90)./(nTrialsYes/2); 
                fa90 = nansum(groupB.(targets{iT}).(fields{iF}).fa90)./(nTrialsYes/2); 
                
                hit90(hit90==1) = 0.99;
                hit90(hit90==0) = 0.01;
                fa90(fa90==1) = 0.99;
                fa90(fa90==0) = 0.01;
                groupB.(targets{iT}).(fields{iF}).sessionDprime90 = norminv(hit90)- norminv(fa90);
                groupB.(targets{iT}).(fields{iF}).subjectDprime90(:,countSubject) = nanmean(norminv(hit90(:,i:i+1)) - norminv(fa90(:,i:i+1)));
            end
            
            countSubject = countSubject + 1;
        end
    end
end

%% plot dprime 

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)    
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iF == iT
            color = p.cueColors(1,:); % valid
        else
            color = p.cueColors(2,:); % invalid
        end
        errorbar(iT,groupB.(targets{iT}).(fields{iF}).dPrime,groupB.(targets{iT}).(fields{iF}).dPrimeSte,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2,'CapSize',12); 
    end
end
xlim([0 3])
xticks([1 2])
xticklabels({'T1' 'T2'})
ylim([0 2.5])
ylabel('sensitivity (d'')') 
if saveFigs
    saveas(gcf,[expt,'_cue_dprime.svg'])
end

%% plot rt

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iT == iF 
            color = p.cueColors(1,:); % valid 
        else
            color = p.cueColors(2,:); % invalid 
        end
        errorbar(iT,nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT,1),2),...
            nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT,1),[],2)/sqrt(numel(subjectNames)),...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
    end
end
xlim([0 3])
xticks([1 2])
xticklabels({'T1' 'T2'})
ylim([0 1.5])
ylabel('reaction time (s)') 
if saveFigs 
    saveas(gcf,[expt,'_cue_rt.svg'])
end

%% plot accuracy 

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iT == iF 
            color = p.cueColors(1,:); % valid 
        else
            color = p.cueColors(2,:); % invalid 
        end
        errorbar(iT,nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc,1),2),...
            nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc,1),[],2)/sqrt(numel(subjectNames)),...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
    end
end
xlim([0 3])
xticks([1 2])
xticklabels({'T1' 'T2'})
ylim([0.5 1])
ylabel('accuracy (proportion correct)') 
if saveFigs 
    saveas(gcf,[expt,'_cue_accuracy.svg'])
end

%% plot dprime x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
                color = p.cueColors(3,:); % neutral
            elseif iT == iF
                color = p.cueColors(1,:); % valid
            else
                color = p.cueColors(2,:); % invalid
            end
            offset = -0.15; 
            mean0 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime0,1),2); 
            ste0 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime0,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean0,ste0,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
        
            offset = 0.15; 
            mean90 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime90,1),2); 
            ste90 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime90,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean90,ste90,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2) 
        end
    end
    xlim([0 3])
    xticks([1 2])
    xticklabels({'T1', 'T2'})
    ylim([0 2.5])
    ylabel('sensitivity (d'')')
    if saveFigs 
        saveas(gcf,[expt,'_cue_dprime_axis.svg'])
    end
end

%% plot rt x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
                color = p.cueColors(3,:); % neutral
            elseif iT == iF
                color = p.cueColors(1,:); % valid
            else
                color = p.cueColors(2,:); % invalid
            end
            offset = -0.15; 
            mean0 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT0,1),2); 
            ste0 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT0,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean0,ste0,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
        
            offset = 0.15; 
            mean90 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT90,1),2); 
            ste90 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectRT90,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean90,ste90,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
        end
    end
    xlim([0 3])
    xticks([1 2])
    xticklabels({'T1', 'T2'})
    ylim([0 1.5])
    ylabel('reaction time (s)')
    if saveFigs
        saveas(gcf,[expt '_cue_rt_axis.svg'])
    end
end

%% plot acc x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
                color = p.cueColors(3,:); % neutral
            elseif iT == iF
                color = p.cueColors(1,:); % valid
            else
                color = p.cueColors(2,:); % invalid
            end
            offset = -0.15; 
            mean0 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc0,1),2); 
            ste0 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc0,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean0,ste0,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
        
            offset = 0.15; 
            mean90 = nanmean(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc90,1),2); 
            ste90 = nanstd(nanmean(groupB.(targets{iT}).(fields{iF}).subjectAcc90,1),[],2)/sqrt(numel(subjectNames)); 
            errorbar(iT+offset,mean90,ste90,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2)
        end
    end
    xlim([0 3])
    xticks([1 2])
    xticklabels({'T1', 'T2'})
    ylim([0.5 1])
    ylabel('accuracy (proportion correct)')
    if saveFigs
        saveas(gcf,[expt '_cue_accuracy_axis.svg'])
    end
end

%% plot dprime sessions

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)    
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iF == iT
            color = p.cueColors(1,:); % valid
        else
            color = p.cueColors(2,:); % invalid
        end
        errorbar(iT,...
            nanmean(groupB.(targets{iT}).(fields{iF}).sessionDprime),...
            nanstd(groupB.(targets{iT}).(fields{iF}).sessionDprime)/sqrt(numel(sessionNames)),...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2,'CapSize',12); 
    end
end
xlim([0 3])
xticks([1 2])
xticklabels({'T1' 'T2'})
ylim([0 2.5])
ylabel('sensitivity (d'')') 
if saveFigs
    saveas(gcf,[expt,'_cue_dprime_sessions.svg'])
end
