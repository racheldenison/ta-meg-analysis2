%% plot behav at subject and group level
% May 2020

%% first run meg_plotBehav for groupB structure 
% input 
targets = {'T1','T2'}; 
expt = 'TANoise'; 
user = 'karen';
slice = 'cueUnsorted'; % cueAcc 'cueAcc_incorrect'
exptDir = meg_pathToTAMEG(expt, user);
% analStr = 'bietfp'; 
saveFigs = 1; 

%% setup
switch expt % conditions 
    case 'TANoise'
        fields = {'cueT1','cueT2'}; 
    case 'TA2'
        fields = {'cueT1','cueT2','neutral'};
end 

[sesssionNames,subjectNames] = meg_sessions(expt); 
subjectsAll = [1 2 3 4 5 6 7 8 9 10]; 
uppers = [1 2 6 8]; 
downers = [3 4 7 9 10]; 

figDir = sprintf('%s/Group/figures/behavior/%s/',exptDir,slice); % figure directory 
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% plot dprime 
figure
hold on 
w = 200; 
h = 500;
pos = [0 0 w h]; % [left bottom width height]
set(gcf,'units','points','position',[100 100 w h]) 
unit = 1.5/numel(subjectNames);
xlim2 = 3; 
for i = 1:numel(subjectNames)*numel(targets) 
    subjectVal = unit*i; 
    xline(subjectVal,'Color',[0.7 0.7 0.7],'LineWidth',0.1)
end

[valSorted,sortIdx] = sort(groupB.T1.cueT1.subjectDprime); % sort by T1 cue neutral dprime 

for iT = 1:numel(targets)
    for iF = 1:numel(fields)  
        vals = groupB.(targets{iT}).(fields{iF}); 
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iF == iT
            color = p.cueColors(1,:); % valid
        else
            color = p.cueColors(2,:); % invalid
        end

        for i = 1:numel(subjectNames)
            if iT == 1
                plot([unit*(i-1),unit*(i)],[vals.subjectDprime(sortIdx(i)),vals.subjectDprime(sortIdx(i))],...
                    'Color',color,'LineWidth',4)
            elseif iT == 2 
                plot([unit*(numel(subjectNames)+i-1),unit*(numel(subjectNames)+i)],...
                    [vals.subjectDprime(sortIdx(i)),vals.subjectDprime(sortIdx(i))],...
                    'Color',color,'LineWidth',4)
            end
        end
        if iT == 1
            plot([0,xlim2/2],[vals.dPrime,vals.dPrime],...
                'Color',color,'LineWidth',4)
            patch([0 0 1.5 1.5],...
                [vals.dPrime-vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime-vals.dPrimeSte],color,'facealpha',0.2,'EdgeColor','none')
        elseif iT ==2
            plot([xlim2/2,xlim2],[vals.dPrime,vals.dPrime],...
                'Color',color,'LineWidth',4) 
            patch([1.5 1.5 3 3],...
                [vals.dPrime-vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime-vals.dPrimeSte],color,'facealpha',0.2,'EdgeColor','none')
        end

%         errorbar(iT,groupB.(targets{iT}).(fields{iF}).dPrime,groupB.(targets{iT}).(fields{iF}).dPrimeSte,...
%             'Color',color,...
%             'Marker','.','MarkerSize',50,...
%             'LineWidth',2,'CapSize',12);
    end
end
 
xlim([0 xlim2])
xticks([xlim2/4 3*xlim2/4])
xticklabels({'T1' 'T2'})
ylim([0 3])
xline(1.5,'k','LineWidth',2)
xline(1.5,'k','LineWidth',2)
ylabel('sensitivity (d'')') 
if saveFigs
    saveas(gcf,sprintf('%s/%s_cue_dprime_individual.svg',figDir,expt))
    saveas(gcf,sprintf('%s/%s_cue_dprime_individual.png',figDir,expt))
end

%% repeat above 
% figure
% hold on 
% w = 200; 
% h = 500;
% pos = [0 0 w h]; % [left bottom width height]
% set(gcf,'units','points','position',[100 100 w h]) 
% unit = 1.5/numel(subjectNames);
% xlim2 = 3; 
% for i = 1:numel(subjectNames)*numel(targets) 
%     subjectVal = unit*i; 
%     xline(subjectVal,'Color',[0.7 0.7 0.7],'LineWidth',0.1)
% end
% 
% [valSorted,sortIdx] = sort(groupB.T1.cueT1.subjectDprime); % sort by T1 cue neutral dprime 
% 
% for iT = 1:numel(targets)
%     for iF = 1:numel(fields)  
%         vals = groupB.(targets{iT}).(fields{iF}); 
%         if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
%             color = p.cueColors(3,:); % neutral
%         elseif iF == iT
%             color = p.cueColors(1,:); % valid
%         else
%             color = p.cueColors(2,:); % invalid
%         end
% 
%         for i = 1:numel(subjectNames)
%             if iT == 1
%                 plot([unit*(i-1),unit*(i)],[vals.subjectDprime(sortIdx(i)),vals.subjectDprime(sortIdx(i))],...
%                     'Color',color,'LineWidth',4)
%             elseif iT == 2 
%                 plot([unit*(numel(subjectNames)+i-1),unit*(numel(subjectNames)+i)],...
%                     [vals.subjectDprime(sortIdx(i)),vals.subjectDprime(sortIdx(i))],...
%                     'Color',color,'LineWidth',4)
%             end
%         end
%         if iT == 1
%             plot([0,xlim2/2],[vals.dPrime,vals.dPrime],...
%                 'Color',color,'LineWidth',4)
%             patch([0 0 1.5 1.5],...
%                 [vals.dPrime-vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime-vals.dPrimeSte],color,'facealpha',0.2,'EdgeColor','none')
%         elseif iT ==2
%             plot([xlim2/2,xlim2],[vals.dPrime,vals.dPrime],...
%                 'Color',color,'LineWidth',4) 
%             patch([1.5 1.5 3 3],...
%                 [vals.dPrime-vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime+vals.dPrimeSte vals.dPrime-vals.dPrimeSte],color,'facealpha',0.2,'EdgeColor','none')
%         end
% 
% %         errorbar(iT,groupB.(targets{iT}).(fields{iF}).dPrime,groupB.(targets{iT}).(fields{iF}).dPrimeSte,...
% %             'Color',color,...
% %             'Marker','.','MarkerSize',50,...
% %             'LineWidth',2,'CapSize',12);
%     end
% end
%  
% xlim([0 xlim2])
% xticks([xlim2/4 3*xlim2/4])
% xticklabels({'T1' 'T2'})
% ylim([0 3])
% xline(1.5,'k','LineWidth',2)
% xline(1.5,'k','LineWidth',2)
% ylabel('sensitivity (d'')') 
% if saveFigs
%     saveas(gcf,sprintf('%s/%s_cue_dprime_individual.svg',figDir,expt))
% end

%% subjects scatter 

figure
hold on
for iF = 1:numel(fields)
    if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
        color = p.cueColors(3,:); % neutral
    elseif iF == iT
        color = p.cueColors(1,:); % valid
    else
        color = p.cueColors(2,:); % invalid
    end
    T1 = groupB.T1.(fields{iF}).subjectDprime;
    T2 = groupB.T2.(fields{iF}).subjectDprime;
    scatter(T1,T2,50,color,'filled')
end
axis equal
axis square
xlim([0 3.5])
ylim([0 3.5])
ylabel('T2 sensitivity (d'')')
xlabel('T1 sensitivity (d'')')
if saveFigs
    saveas(gcf,sprintf('%s/%s_cue_dprime_scatter.svg',figDir,expt))
    saveas(gcf,sprintf('%s/%s_cue_dprime_scatter.png',figDir,expt))
end

%% dprime difference sorted 
figure
hold on 
w = 200; 
h = 500;
pos = [0 0 w h]; % [left bottom width height]
set(gcf,'units','points','position',[100 100 w h]) 
unit = 1.5/numel(subjectNames);
xlim2 = 3; 
for i = 1:numel(subjectNames)*numel(targets) 
    subjectVal = unit*i; 
    xline(subjectVal,'Color',[0.7 0.7 0.7],'LineWidth',0.1)
end

[valSorted,sortIdx] = sort(groupB.T1.cueT1.subjectDprime-groupB.T1.neutral.subjectDprime,'descend'); % sort by T1 cue neutral dprime 

for iT = 1:numel(targets)
    for iF = 1:2 % valid, invalid, no neutral 
        vals = groupB.(targets{iT}).(fields{iF}); 
        neutral = groupB.(targets{iT}).neutral; 
        if iF == iT
            color = p.cueColors(1,:); % valid
            change = vals.subjectDprime - neutral.subjectDprime; 
        else
            color = p.cueColors(2,:); % invalid
            change = -(neutral.subjectDprime - vals.subjectDprime); 
        end

        for i = 1:numel(subjectNames)
            if iT == 1
                plot([unit*(i-1),unit*(i)],[change(sortIdx(i)),change(sortIdx(i))],...
                    'Color',color,'LineWidth',4)
            elseif iT == 2 
                plot([unit*(numel(subjectNames)+i-1),unit*(numel(subjectNames)+i)],...
                    [change(sortIdx(i)),change(sortIdx(i))],...
                    'Color',color,'LineWidth',4)
            end
        end
        hline(0,':k')
        meanChange = nanmean(change); 
        steChange = nanstd(change)/sqrt(numel(subjectNames)); 
        if iT == 1
            plot([0,xlim2/2],[meanChange,meanChange],...
                'Color',color,'LineWidth',4)
            patch([0 0 1.5 1.5],...
                [meanChange-steChange meanChange+steChange meanChange+steChange meanChange-steChange],color,'facealpha',0.2,'EdgeColor','none')
        elseif iT ==2
            plot([xlim2/2,xlim2],[meanChange,meanChange],...
                'Color',color,'LineWidth',4) 
            patch([1.5 1.5 3 3],...
                [meanChange-steChange meanChange+steChange meanChange+steChange meanChange-steChange],color,'facealpha',0.2,'EdgeColor','none')
        end
    end
end
 
xlim([0 xlim2])
xticks([xlim2/4 3*xlim2/4])
xticklabels({'T1' 'T2'})
ylim([-1.5 1.5])
xline(1.5,'k','LineWidth',2)
xline(1.5,'k','LineWidth',2)
ylabel('Δ sensitivity from neutral (d'')') 
set(gca,'TickDir','out');
if saveFigs
    % saveas(gcf,[expt,'_cue_delta_dprime_individual.svg'])
end

%% scatterplot T1 valid, T2 invalid 

figure
hold on 
w = 200; 
h = 500;
pos = [0 0 w h]; % [left bottom width height]
set(gcf,'units','points','position',[100 100 w h]) 
unit = 1.5/numel(subjectNames);
xlim2 = 3; 
for i = 1:numel(subjectNames)*numel(targets) 
    subjectVal = unit*i; 
    xline(subjectVal,'Color',[0.7 0.7 0.7],'LineWidth',0.1)
end

[valSorted,sortIdx] = sort(groupB.T1.cueT1.subjectDprime-groupB.T1.neutral.subjectDprime,'descend'); % sort by T1 cue neutral dprime 

scatter(groupB.T1.cueT1.subjectDprime,groupB.T2.cueT1.subjectDprime)
hold on 
scatter(groupB.T2.cueT2.subjectDprime,groupB.T1.cueT2.subjectDprime)
hold on 
scatter(groupB.T1.cueT1.subjectDprime,groupB.T1.cueT2.subjectDprime)

xlim([0 xlim2])
xticks([xlim2/4 3*xlim2/4])
xticklabels({'T1' 'T2'})
ylim([-1.5 1.5])
xline(1.5,'k','LineWidth',2)
xline(1.5,'k','LineWidth',2)
ylabel('Δ sensitivity from neutral (d'')') 
set(gca,'TickDir','out');
if saveFigs
    % saveas(gcf,[expt,'_cue_dprime.svg'])
end

%% plot individual __ variable 

sorted = 0; 
variable = 'subjectDprime'; % 'subjectRT'; % 'subjectDprime', 'subjectAcc'

figure
hold on 
w = 200; 
h = 500;
pos = [0 0 w h]; % [left bottom width height]
set(gcf,'units','points','position',[100 100 w h]) 
unit = 1.5/numel(subjectNames);
xlim2 = 3; 
for i = 1:numel(subjectNames)*numel(targets) 
    subjectVal = unit*i; 
    xline(subjectVal,'Color',[0.7 0.7 0.7],'LineWidth',0.1)
end

if sorted
    [valSorted,sortIdx] = sort(groupB.T1.cueT1.(variable)); % sort by T1 cue neutral dprime 
end

for iT = 1:numel(targets)
    for iF = 1:numel(fields)  
        vals = groupB.(targets{iT}).(fields{iF}); 
        if strcmp(expt,'TA2') && strcmp(fields{iF},'neutral')
            color = p.cueColors(3,:); % neutral
        elseif iF == iT
            color = p.cueColors(1,:); % valid
        else
            color = p.cueColors(2,:); % invalid
        end

        for i = 1:numel(subjectNames)
            if iT == 1
                if sorted
                    plot([unit*(i-1),unit*(i)],[vals.(variable)(sortIdx(i)),vals.(variable)(sortIdx(i))],...
                        'Color',color,'LineWidth',4)
                else
                    plot([unit*(i-1),unit*(i)],[vals.(variable)(i),vals.(variable)(i)],...
                        'Color',color,'LineWidth',4)
                end
            elseif iT == 2
                if sorted
                    plot([unit*(numel(subjectNames)+i-1),unit*(numel(subjectNames)+i)],...
                        [vals.(variable)(sortIdx(i)),vals.(variable)(sortIdx(i))],...
                        'Color',color,'LineWidth',4)
                else
                    plot([unit*(numel(subjectNames)+i-1),unit*(numel(subjectNames)+i)],...
                        [vals.(variable)(i),vals.(variable)(i)],...
                        'Color',color,'LineWidth',4)
                end
            end
        end
        varMean = nanmean(vals.(variable));
        varSte = nanstd(vals.(variable))/sqrt(numel(subjectNames)); 
        if iT == 1
            plot([0,xlim2/2],[varMean,varMean],...
                'Color',color,'LineWidth',4)
            patch([0 0 1.5 1.5],...
                [varMean-varSte varMean+varSte varMean+varSte varMean-varSte],color,'facealpha',0.2,'EdgeColor','none')
        elseif iT ==2
            plot([xlim2/2,xlim2],[varMean,varMean],...
                'Color',color,'LineWidth',4) 
            patch([1.5 1.5 3 3],...
                [varMean-varSte varMean+varSte varMean+varSte varMean-varSte],color,'facealpha',0.2,'EdgeColor','none')
        end
    end
end
 
xlim([0 xlim2])
xticks([xlim2/4 3*xlim2/4])
xticklabels({'T1' 'T2'})
% ylim([0 3])
xline(1.5,'k','LineWidth',2)
xline(1.5,'k','LineWidth',2)
switch variable 
    case 'subjectRT'
        ylabel('reaction time (ms)') 
    case 'subjectDprime'
        ylabel('sensitivity (d'')') 
    case 'subjectAcc'
        ylabel('accuracy (proportion correct)') 
end
if saveFigs
    saveas(gcf,sprintf('%s/%s_cue_%s.svg',figDir,expt,variable))
    saveas(gcf,sprintf('%s/%s_cue_%s.png',figDir,expt,variable))
end




