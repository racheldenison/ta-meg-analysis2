function meg_plotBehav(expt, user)

% function meg_plotBehav(expt, user)
% plots sensitivity, accuracy, and rt by condition by axis orientation
% (optional) 
%
% INPUTS
% expt
%   string, 'TA2' or 'TANoise'
% user (optional)
%   string giving the user name to determine path to data. defaults to
%   'mcq' = get the data from the mcq server
%
% Karen Tian
% May 2020

%% inputs
expt = 'TA2'; % TANoise 
user = 'karen'; 

[sessionNames,subjectNames] = meg_sessions(expt); 
targets = {'T1','T2'}; 
plotAxis = 1; % slice by axis orientation 

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

%% collect behavior all sessions 
groupB = [];
for i = 1:numel(sessionNames)
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            targetIdx = find(groupDAll(i).B.responseTarget==iT);
            fieldIdx = groupDAll(i).I.(fields{iF}); 
            unionIdx = intersect(targetIdx,fieldIdx); 
            
            groupB.(targets{iT}).(fields{iF}).acc(:,i) = groupDAll(i).B.acc(unionIdx);
            groupB.(targets{iT}).(fields{iF}).rt(:,i) = groupDAll(i).B.rt(unionIdx);
            
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

%% avg behavior subject level 

for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        nTrials = size(groupB.(targets{iT}).(fields{iF}).acc,1); 
        countSubject = 1;
        for i = 1:2:numel(sessionNames)
            groupB.(targets{iT}).(fields{iF}).subjectAcc(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).acc(:,i:i+1),2);
            groupB.(targets{iT}).(fields{iF}).subjectRT(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).rt(:,i:i+1),2);
            
            groupB.(targets{iT}).(fields{iF}).subjectHit(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).hit(:,i:i+1),2);
            groupB.(targets{iT}).(fields{iF}).subjectFA(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).fa(:,i:i+1),2);
            
            % dprime 
            hit = nansum(groupB.(targets{iT}).(fields{iF}).subjectHit(:,countSubject))/nTrials; 
            fa = nansum(groupB.(targets{iT}).(fields{iF}).subjectFA(:,countSubject))/nTrials; 
            
            hit(hit==1) = 0.99; 
            hit(hit==0) = 0.01; 
            fa(fa==1) = 0.99; 
            fa(fa==0) = 0.01;
            
            groupB.(targets{iT}).(fields{iF}).subjectDprime(countSubject) = norminv(hit)- norminv(fa);
            
            if plotAxis 
                groupB.(targets{iT}).(fields{iF}).subjectAcc0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).acc0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectAcc90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).acc90(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectRT0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).rt0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectRT90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).rt90(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectHit0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).hit0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectHit90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).hit90(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectFA0(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).fa0(:,i:i+1),2);
                groupB.(targets{iT}).(fields{iF}).subjectFA90(:,countSubject) = nanmean(groupB.(targets{iT}).(fields{iF}).fa90(:,i:i+1),2);
                
                hit0 = nansum(groupB.(targets{iT}).(fields{iF}).subjectHit0(:,countSubject))/(nTrials/2);
                fa0 = nansum(groupB.(targets{iT}).(fields{iF}).subjectFA0(:,countSubject))/(nTrials/2);
                hit0(hit0==1) = 0.99;
                hit0(hit0==0) = 0.01;
                fa0(fa0==1) = 0.99;
                fa0(fa0==0) = 0.01;
                groupB.(targets{iT}).(fields{iF}).subjectDprime0(countSubject) = norminv(hit0)- norminv(fa0);
                
                hit90 = nansum(groupB.(targets{iT}).(fields{iF}).subjectHit90(:,countSubject))/(nTrials/2);
                fa90 = nansum(groupB.(targets{iT}).(fields{iF}).subjectFA90(:,countSubject))/(nTrials/2);
                hit90(hit90==1) = 0.99;
                hit90(hit90==0) = 0.01;
                fa90(fa90==1) = 0.99;
                fa90(fa90==0) = 0.01;
                groupB.(targets{iT}).(fields{iF}).subjectDprime90(countSubject) = norminv(hit90)- norminv(fa90);
            end
            
            countSubject = countSubject + 1;
        end
        groupB.(targets{iT}).(fields{iF}).dPrime = nanmean(groupB.(targets{iT}).(fields{iF}).subjectDprime);
        groupB.(targets{iT}).(fields{iF}).dPrimeSte = nanstd(groupB.(targets{iT}).(fields{iF}).subjectDprime)/sqrt(numel(subjectNames));
    end
end

%% plot dprime 

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        if iT == iF 
            color = p.cueColors(1,:); % valid 
        else
            color = p.cueColors(2,:); % invalid 
        end
        errorbar(iT,groupB.(targets{iT}).(fields{iF}).dPrime,groupB.(targets{iT}).(fields{iF}).dPrimeSte,...
            'Color',color,...
            'Marker','.','MarkerSize',30,...
            'LineWidth',2); 
    end
end
xlim([0 3])
xticks([1 2])
xticklabels({'T1' 'T2'})
ylabel('sensitivity (d'')') 
% saveas(gcf,'cue_dprime.svg')

%% plot rt

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        if iT == iF 
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
ylabel('reaction time (s)') 
% saveas(gcf,'cue_rt.svg')

%% plot proportion correct 

figure
hold on 
set(gcf,'Position',[100 100 200 400]) 
for iT = 1:numel(targets)
    for iF = 1:numel(fields)
        if iT == iF 
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
ylabel('accuracy (proportion correct)') 
% set(gca,'fontname','helvetica','fontsize',14)
% saveas(gcf,'cue_accuracy.svg')

%% plot dprime x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if iT == iF
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
    ylabel('sensitivity (d'')')
    saveas(gcf,'cue_dprime_axis.svg')
end

%% plot rt x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if iT == iF
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
    ylabel('reaction time (s)')
    saveas(gcf,'cue_rt_axis.svg')
end

%% plot proportion x axis 

if plotAxis
    figure
    hold on
    set(gcf,'Position',[100 100 200 400])
    for iT = 1:numel(targets)
        for iF = 1:numel(fields)
            if iT == iF
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
    ylabel('accuracy (proportion correct)')
    % saveas(gcf,'cue_accuracy_axis.svg')
end

