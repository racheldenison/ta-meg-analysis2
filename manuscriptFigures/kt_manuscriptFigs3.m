% ITPC spectrogram by cue 
% Average trial power 

% settings 
p = meg_params('TANoise_ITPCsession8'); 
saveFigs = 0; 
figDir = '/Users/kantian/Dropbox/Data/TANoise/Manuscript_Figs'; 

%% figure ITPC spectrogram by cue 
load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/ITPC_spectrogram/groupA_ITPCspectrogram_byAtt.mat')

toi = -2000:5000; 
foi = 1:50; 
ytick = 10:10:numel(foi);
xtick = -2000:500:5000;

ylims = [min(foi),max(foi)]; 
xlims = [abs(p.tstart)-100,abs(p.tstart)+2400]; 

% cue T1 
val = mean(A.cueT1.subject,3); 
figure
set(gcf,'Position',[100 100 500 400])
imagesc(val) 
meg_timeFreqPlotLabels2(toi,foi,xtick,ytick,p.eventTimes+abs(p.tstart))
xlim(xlims)
ylim(ylims)
meg_figureStyle
title('Cue T1') 
xticks(2000:500:4000)
xticklabels({'0','500','1000','1500','2000'})
c = colorbar;
caxis([0.05 0.35])
xlabel('Time (ms)') 
ylabel('Frequency') 
c.Box = 'off'; 
c.Label.String = 'ITPC (A.U.)'; 
if saveFigs
    figTitle = 'ITPCSpectrogram_CueT1'; 
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

% cue T2
val = mean(A.cueT2.subject,3); 
figure
set(gcf,'Position',[100 100 500 400])
imagesc(val) 
meg_timeFreqPlotLabels2(toi,foi,xtick,ytick,p.eventTimes+abs(p.tstart))
xlim(xlims)
ylim(ylims)
meg_figureStyle
title('Cue T2') 
xticks(2000:500:4000)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Time (ms)') 
ylabel('Frequency') 
c = colorbar;
caxis([0.05 0.35])
c.Box = 'off'; 
c.Label.String = 'ITPC (A.U.)'; 
if saveFigs
    figTitle = 'ITPCSpectrogram_CueT2'; 
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

% cue T1 - cue T2
val = mean(A.cueT1.subject-A.cueT2.subject,3); 
figure
set(gcf,'Position',[100 100 500 400])
imagesc(val) 
meg_timeFreqPlotLabels2(toi,foi,xtick,ytick,p.eventTimes+abs(p.tstart))
xlim(xlims)
ylim(ylims)
meg_figureStyle
title('Cue T1 - Cue T2') 
xticks(2000:500:4000)
xticklabels({'0','500','1000','1500','2000'})
xlabel('Time (ms)') 
ylabel('Frequency') 
c = colorbar;
caxis([-0.045 0.045])
c.Ticks = linspace(-0.045,0.045,7);
c.Box = 'off'; 
c.Label.String = 'ITPC (A.U.)'; 
if saveFigs
    figTitle = 'ITPCSpectrogram_CueT1-T2'; 
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle)) 
end

%% Average trial power - normalize 
% average all trials data session to subjects --> baseline -->  save
for i = 1:20
    val = []; 
    val = groupTFspectrogram_avgTrial(i).pows(20,:); 
    vals(i,:) = val; 
end
% average to subject 
count = 1; 
for i = 1:2:19
    val = []; 
    val = vals(i:i+1,:); 
    val = mean(val,1); 
    vals_sub(count,:) = val; 
    count = count + 1; 
end

% baseline mat 
normalizeToi = 47:97; 
normalizeIdx = abs(p.tstart)/100 + normalizeToi; 
baselineMat = vals_sub(:,normalizeIdx); 
baselineMat = mean(baselineMat,2); 

% subtract baseline mat and get itpc peak direction 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions('TANoise'); 
for i = 1:10
    normPow(i,:) = (vals_sub(i,:) - baselineMat(i)) * ITPCsubject(i); 
end 

% average to group 
normPow_group = mean(normPow,1,'omitnan'); 

% plot norm power group 
figure 
plot(normPow_group)

%% Power by cue 
vals_cueT1 = squeeze(A.cueT1.subject(20,:,:)); % subjects 
vals_cueT2 = squeeze(A.cueT2.subject(20,:,:)); 

% subtract baseline mat and flip 
for i = 1:10
    normVals_cueT1(:,i) = (vals_cueT1(:,i) - baselineMat(i)) * ITPCsubject(i); 
    normVals_cueT2(:,i) = (vals_cueT2(:,i) - baselineMat(i)) * ITPCsubject(i); 
end

g_normVals_cueT1 = mean(normVals_cueT1,2); 
g_normVals_cueT2 = mean(normVals_cueT2,2); 
% plot 
figure
hold on 
plot(g_normVals_cueT1,'LineWidth',2)
plot(g_normVals_cueT2,'LineWidth',2)
%% 
t = -2000:10:5000; 

figure
hold on 
plot(t,A.cueT1.group(20,:)) 
plot(t,A.cueT2.group(20,:)) 
meg_figureStyle
xlim([-100 2500])










