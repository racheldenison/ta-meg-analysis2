%% generate manuscript figures
% April 2021

%% setup 
expt = 'TANoise'; 
p = meg_params('TANoise_ITPCsession8'); 

[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 
t = p.tstart:p.tstop; 
toi = p.eventTimes(1):p.eventTimes(4); % time of interest: precue -> response cue 
tBuffer = 100; 
toiBuffer = toi(1)-tBuffer:toi(end)+tBuffer; 
xlims = [size(toiBuffer,1),size(toiBuffer,2)]; 
% xlims = [1,length(toiBuffer)]; 
foi = 1:50; 
ylims = [foi(1),foi(end)]; 

ytick = 10:10:numel(foi);
xtick = tBuffer+1:500:length(toi)+tBuffer+1;

%% 
saveFigs = 1; 
figDir = sprintf('%s/figures',pwd); 
disp(figDir)
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% load data %%%
% A
% load('/Users/kantian/Dropbox/Data/TANoise/MEG/Group/mat/ITPC_spectrogram/groupA_ITPCspectrogram_byAtt.mat'); 
load('groupA_ITPCspectrogram_byAtt.mat')

%% ITPC spectrogram (cueT1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable = 'normSubject'; % 'normSubjectFlipped' 'normSubject'
% vals 
vals_cueT1 = A.cueT1.(variable)(:,find(t==toiBuffer(1)):find(t==toiBuffer(end)),:); 
vals_cueT2 = A.cueT2.(variable)(:,find(t==toiBuffer(1)):find(t==toiBuffer(end)),:); 
% 

vals = vals_cueT1; 
vals = mean(vals,3,'omitnan'); % average across subjects 

figure
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
switch variable 
    case 'normSubjectFlipped'
        caxis([-0.1 0.1])
    case 'normSubject'
        caxis([-0.1 0.2])
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Cue T1'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_cueT1.svg', figDir, variable)) 
end

%% ITPC spectrogram (cueT2)  
vals = vals_cueT2; 
vals = mean(vals,3,'omitnan'); % average across subjects 

figure
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
switch variable 
    case 'normSubjectFlipped'
        caxis([-0.1 0.1])
    case 'normSubject'
        caxis([-0.1 0.2])
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Cue T2'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_cueT2.svg', figDir, variable)) 
end

%% ITPC spectrogram (cueT1 - cueT2) 
vals = vals_cueT1 - vals_cueT2; 
vals = mean(vals,3,'omitnan'); % average across subjects 

% permutation test 

figure
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
caxis([-0.04 0.04])
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Cue T1 - Cue T2'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_cueT1-cueT2.svg', figDir, variable)) 
end

%% ITPC spectrogram (uppers - downers) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable = 'normSubject'; % 'normSubject' 'normSubjectFlipped' 'subject'

idxUppers = find(ITPCsubject==1); 
idxDowners = find(ITPCsubject==-1); 

vals_uppers = A.all.(variable)(:,find(t==toiBuffer(1)):find(t==toiBuffer(end)),idxUppers); 
vals_downers = A.all.(variable)(:,find(t==toiBuffer(1)):find(t==toiBuffer(end)),idxDowners); 

%% ITPC spectrogram uppers 
vals = vals_uppers; 
vals = mean(vals,3,'omitnan'); 

figure 
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
switch variable 
    case 'normSubjectFlipped'
        % caxis([-0.1 0.1])
    case 'normSubject'
        caxis([-0.1 0.2])
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Uppers'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_uppers.svg', figDir, variable)) 
end

%% ITPC spectrogram downers
vals = vals_downers; 
vals = mean(vals,3,'omitnan'); 

figure 
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
switch variable 
    case 'normSubjectFlipped'
        % caxis([-0.1 0.1])
    case 'normSubject'
        caxis([-0.1 0.2])
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Downers'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_downers.svg', figDir, variable)) 
end

%% ITPC spectrogram (uppers - downers)
vals = mean(vals_uppers,3,'omitnan') - mean(vals_downers,3,'omitnan'); 

figure 
set(gcf, 'Position',  [100, 100, 500, 500])
imagesc(vals) 
xlim(xlims)
ylim(ylims)
switch variable 
    case 'normSubjectFlipped'
        % caxis([-0.4 0.4])
    case 'normSubject'
        caxis([-0.15 0.15])
end
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
cbar = colorbar; 
cbar.Label.String = 'ITPC'; 
meg_timeFreqPlotLabels(toiBuffer,foi,xtick,ytick,p.eventTimes+tBuffer);
meg_figureStyle
title(sprintf('Uppers - Downers'))

if saveFigs 
    saveas(gcf,sprintf('%s/ITPCspectrogram_%s_uppers-downers.svg', figDir, variable)) 
end

