% December 21, 2021 

%% save baseline ITPC 
foi = 20; 
baselineVals = A.all.baselineITPC(foi,:); % foi x sessions 

% average to subject 
count = 1; 
for i = 1:2:19
    vals = []; 
    vals = baselineVals(i:i+1); 
    baselineVals_subject(count) = mean(vals); 
    count = count+1; 
end

%% plot peak T1 ITPC 
% --- Figure ---
figure
set(gcf,'Position',[100 100 300 400])
hold on 
meg_figureStyle
xJitter = 0.2; 
idxPeaks = 1:20; 
idxPeaks(10:11) = []; 

peak_T1 = A.cueT1.peakT1ITPC(idxPeaks);
peak_T2 = A.cueT2.peakT1ITPC(idxPeaks);

% T1 
% valid 
x_valid = repmat(1-xJitter,[18,1]);
scatter(x_valid,peak_T1,'filled','MarkerFaceColor',p.cueColors(1,:),'MarkerFaceAlpha',0.5)
errorbar(x_valid(1),mean(peak_T1),std(peak_T1)/sqrt(20),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(1,:),'MarkerEdgeColor',p.cueColors(1,:),...
    'Color',p.cueColors(1,:),'LineWidth',2);
% invalid 
x_valid = repmat(1+xJitter,[18,1]);
scatter(x_valid,peak_T2,'filled','MarkerFaceColor',p.cueColors(2,:),'MarkerFaceAlpha',0.5)
errorbar(x_valid(1),mean(peak_T2),std(peak_T2)/sqrt(20),'Marker','.','MarkerSize',40,'MarkerFaceColor',p.cueColors(2,:),'MarkerEdgeColor',p.cueColors(2,:),...
    'Color',p.cueColors(2,:),'LineWidth',2); 

xlim([0 3])