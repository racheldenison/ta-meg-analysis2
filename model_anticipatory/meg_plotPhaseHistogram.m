function [group,fittedP,stats] = meg_plotPhaseHistogram(mdlFit,paramNames,figDir,mdlFitType)
% function meg_plotPhaseHistogram(mdlFit,figDir,mdlFitType)
% mdlFitType can be 'separate' 'yoked_fixedFreq' 'yoked'

%% Setup
saveFigs = 1;
edges = 0:pi/24:2*pi; % 24 bins, try 12 bins -pi to pi? to match the bounds? 
fitLevel = 'session'; 

% -- MEG settings --- 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 
Fs = 1000; % sampling frequency per second 
% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% -- Export CSV --- 
exportCSV = 0; 
csvDir = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data';
if ~exist(csvDir, 'dir')
    mkdir(csvDir)
end

%% --- Model parameters --- 
clear fvals solutions idxPhase 
switch mdlFitType
    case 'separate' % fixed frequency
        idxPhase = [find(contains(paramNames,'phase'))];
        % --- Precue T1 ---
        fvals(:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2) 
        solutions(:,:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2) 
        % --- Precue T2 ---
        fvals(:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2) 
        solutions(:,:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2) 
        
        idxAmp = [find(contains(paramNames,'amplitude'))];
        
        freq = mdlFit.linear2Hz.cueT1.(fitLevel).fixedFreq;  
        fprintf('Frequency fixed %d Hz',freq)
    case 'yoked_fixedFreq'
        idxPhase = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))];
        fvals = mdlFit.linear2Hz.(fitLevel).fval;
        solutions = mdlFit.linear2Hz.(fitLevel).solution; % permutations (100) x subjects (20) x parameters (8) 
        freq = 2; 
        fprintf('Frequency fixed %d Hz',freq)
    case 'yoked'
        idxPhase = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))];
        idxFreq = [find(contains(paramNames,'freq'))]; 
        fvals = mdlFit.linear2Hz.(fitLevel).fval;
        solutions = mdlFit.linear2Hz.(fitLevel).solution; % permutations (100) x subjects (20) x parameters (9) 
end
nPerms = size(fvals,1); 
nSubs = size(fvals,2); 
nParams = size(solutions,3); 
rhoPerm = 30; 
cueLevel = {'cueT1','cueT2'};

%% --- Extract amplitude 
clear valAmp1_min valAmp2_min
for iS = 1:nSubs
    switch mdlFitType
        case 'separate'
            % --- Precue T1 ---
            [minVal1,idx1] = min( fvals(:,iS,1) );
            valAmp1 = solutions(:,iS,idxAmp,1); % amplitude  
            valAmp1_min(iS) = valAmp1(idx1);
            % --- Precue T2 ---
            [minVal2,idx2] = min( fvals(:,iS,2) );
            valAmp2 = solutions(:,iS,idxAmp,2); % amplitude  
            valAmp2_min(iS) = valAmp2(idx2);
        case 'yoked'
            [minVal,idx] = min( fvals(:,iS) );
            idxAmp1 = [find(contains(paramNames,'amplitude1'))];
            idxAmp2 = [find(contains(paramNames,'amplitude2'))];
            valAmp1 = solutions(:,iS,idxAmp1); % amplitude precue T1 
            valAmp1_min(iS) = valAmp1(idx);
            valAmp2 = solutions(:,iS,idxAmp2); % amplitude precue T2 
            valAmp2_min(iS) = valAmp2(idx);
    end
end
amp.valAmp1_min = valAmp1_min; 
amp.valAmp2_min = valAmp2_min; 

%% --- Extract frequency (for yoked free freq) --- 
switch mdlFitType
    case 'yoked'
        idxFreq = [find(contains(paramNames,'freq'))];
        for iS = 1:nSubs
            [minVal,idx] = min( fvals(:,iS) );
            valFreq = solutions(:,iS,idxFreq); % frequency 
            valFreq_min(iS) = valFreq(idx);
        end
end

% -- Histogram of fitted freqs (sessions) --- 

%% --- Extract and plot fitted phases (sessions) --- 
for iS = 1:nSubs
    figure
    set(gcf,'Position',[100 100 800 200])
    switch mdlFitType
        case 'separate'
            % --- Precue T1 ---
            [minVal1,idx1] = min( fvals(:,iS,1) );
            valR1 = solutions(:,iS,idxPhase,1); % radian space 
            valR1_min = valR1(idx1);
            valT1_min = rad2t(valR1_min,freq,Fs); % t space (s/100) 
            % --- Precue T2 ---
            [minVal2,idx2] = min( fvals(:,iS,2) );
            valR2 = solutions(:,iS,idxPhase,2); % radian space 
            valR2_min = valR2(idx2); 
            valT2_min = rad2t(valR2_min,freq,Fs); 
            % --- Phase difference ---
            valR_diffAbs = abs(circ_dist(valR1,valR2));
            valR_diffMinAbs = abs(circ_dist(valR1_min,valR2_min)); 
            valT_diffMinAbs = rad2t(valR_diffMinAbs,freq,Fs); 
        case 'yoked_fixedFreq'
            [minVal,idx] = min( fvals(:,iS) );
             % --- Precue T1 ---
            valR1 = solutions(:,iS,idxPhase(1)); % radian space 
            valR1_min = valR1(idx);
            valT1_min = rad2t(valR1_min,freq,Fs); % t space (s/100) 
            % --- Precue T2 ---
            valR2 = solutions(:,iS,idxPhase(2)); % radian space 
            valR2_min = valR2(idx); 
            valT2_min = rad2t(valR2_min,freq,Fs);
        case 'yoked'
            [minVal,idx] = min( fvals(:,iS) );
            freq = solutions(idx,iS,idxFreq); % Hz
            % --- Precue T1 ---
            valR1 = solutions(:,iS,idxPhase(1)); % radian space
            valR1_min = valR1(idx);
            valT1_min = rad2t(valR1_min,freq,Fs); % t space (s/100)
            % --- Precue T2 ---
            valR2 = solutions(:,iS,idxPhase(2)); % radian space
            valR2_min = valR2(idx);
            valT2_min = rad2t(valR2_min,freq,Fs);
    end
    % --- Phase difference ---
    valR_diff = circ_dist(valR1,valR2); 
    valR_diffAbs = abs(circ_dist(valR1,valR2));

    valR_diffMin = circ_dist(valR1_min,valR2_min); 
    valT_diffMin = rad2t(valR_diffMin,freq,Fs); 
    valR_diffMinAbs = abs(circ_dist(valR1_min,valR2_min));
    valT_diffMinAbs = rad2t(valR_diffMinAbs,freq,Fs); 

    % --- Plot precue T1 --- 
    subplot 141
    polarhistogram(valR1,edges); % fitted phase histogram
    hold on 
    polarplot([0 valR1_min],[0 rhoPerm],'r') % fitted phase with min fval 
    txt = sprintf('Precue T1\n%0.2f ms/100 = %0.2f\\pi rad',...
        valT1_min,...
        valR1_min/pi);
    title(txt)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';

    % --- Plot precue T2 --- 
    subplot 142
    polarhistogram(valR2,edges); % fitted phase histogram
    hold on 
    polarplot([0 valR2_min],[0 rhoPerm],'r') % fitted phase with min fval 
    txt = sprintf('Precue T2\n%0.2f ms/100 = %0.2f\\pi rad',...
        valT2_min,...
        valR2_min/pi);
    title(txt)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';

    % --- Plot precue T1 - precue T2 difference ---
    subplot 143
    polarhistogram(valR_diff,edges); % fitted phase histogram
    hold on
    polarplot([0 valR_diffMin],[0 rhoPerm],'r') % fitted phase with min fval
    txt = sprintf('Precue T1 - Precue T2\n%0.2f ms/100 = %0.2f\\pi rad',...
        valT_diffMin,...
        valR_diffMin/pi);
    title(txt)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';

    % --- Plot precue T1 - precue T2 abs difference ---
    subplot 144
    polarhistogram(valR_diffAbs,edges); % fitted phase histogram
    hold on
    polarplot([0 valR_diffMinAbs],[0 rhoPerm],'r') % fitted phase with min fval
    txt = sprintf('abs(Precue T1 - Precue) T2\n%0.2f ms/100 = %0.2f\\pi rad',...
        valT_diffMinAbs,...
        valR_diffMinAbs/pi);
    title(txt)
    pax = gca;
    pax.ThetaAxisUnits = 'radians';

    sgtitle(sprintf('%s, nPerms = %d',...
        und2space(sessionNames{iS}),...
        nPerms))

    if saveFigs
        figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_PermutationsHistogram_%d',...
            sessionNames{iS},...
            mdlFitType,...
            nPerms);
        saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    end
    
    % save to group structure 
    group.valR1_min(iS) = valR1_min;
    group.valT1_min(iS) = valT1_min;

    group.valR2_min(iS) = valR2_min;
    group.valT2_min(iS) = valT2_min;

    group.valR_diffMin(iS) = valR_diffMin;
    group.valT_diffMin(iS) = valT_diffMin;

    group.valR_diffMinAbs(iS) = valR_diffMinAbs;
    group.valT_diffMinAbs(iS) = valT_diffMinAbs;
end

%% Plot group (phase) --- session level --- 
figure
set(gcf,'Position',[100 100 800 200])

% --- Plot precue T1 ---
subplot 141
polarhistogram(group.valR1_min,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T1')

% --- Plot precue T2 ---
subplot 142
polarhistogram(group.valR2_min,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T2')

% --- Plot precue T1 - precue T2 difference ---
subplot 143
polarhistogram(group.valR_diffMin,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T1 - T2')

% --- Plot precue T1 - precue T2 abs difference ---
subplot 144
polarhistogram(group.valR_diffMinAbs,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('abs(Precue T1 - T2)')

sgtitle(sprintf('Group (n = %d)',...
    nSubs))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_PermutationsHistogram_%d',...
        mdlFitType,...
        nSubs);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Plot group (phase) by amplidue --- Sessions --- 
figure
set(gcf,'Position',[100 100 800 200])

% --- Plot precue T1 ---
subplot 141
% rho is amplitude 
rho = valAmp1_min; 
% theta is phase 
theta = group.valR1_min; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T1')

% --- Plot precue T2 ---
subplot 142
% rho is amplitude 
rho = valAmp2_min; 
% theta is phase 
theta = group.valR2_min; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T2')

% --- Plot precue T1 - precue T2 difference ---
subplot 143
% rho is amplitude difference of precue T1 and T2 
rho = valAmp1_min-valAmp2_min; 
% theta is phase 
theta = group.valR_diffMin; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1 - T2\nAmp difference'))

% --- Plot precue T1 - precue T2 difference ---
subplot 144
% rho is amplitude mean of precue T1 and T2 
rho = (valAmp1_min+valAmp2_min)/2; 
% theta is phase 
theta = group.valR_diffMin; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1 - T2\nAmp mean'))

sgtitle(sprintf('Group (n = %d) Phase x Amplitude',...
    nSubs))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear%dHz_Phase-Amp_Sessions%d',...
        mdlFitType,freq,...
        nSubs);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Plot group (phase) --- session level (OVERLAID) --- 
figure
set(gcf,'Position',[100 100 300 200])

val = group.valR1_min; % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(1,:),'FaceAlpha',1); % fitted phase histogram
hold on
val = group.valR2_min; % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(2,:),'FaceAlpha',0.7); 
pax = gca;
pax.ThetaAxisUnits = 'radians';

sgtitle(sprintf('Sessions (n = 20)'))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_PermutationsHistogram_sessions_overlaid',...
        mdlFitType);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Compile phase difference from group struct 
clear fittedP 
subCount = 0;
fields = fieldnames(group); 
for iS = 1:2:20 % sessions
    subCount = subCount+1;
    for iV = 1:numel(fields)
        clear val
        val = group.(fields{iV})(iS:iS+1); 
        fittedP.linear2Hz.(fields{iV})(subCount,:) = val; % subjects(10) x sessions (2) 
    end
end
% group.subjects = fittedP; 

% Add amplitude 
ampFields = {'valAmp1_min','valAmp2_min'};
subCount = 0;
for iS = 1:2:20 % sessions
    subCount = subCount+1;
    for iV = 1:numel(ampFields)
        clear val
        val = amp.(ampFields{iV})(iS:iS+1);
        fittedP.linear2Hz.(ampFields{iV})(subCount,:) = val;
    end
end

%% Plot group (freq) 
switch mdlFitType
    case 'yoked'
end

%% Plot group (phase) --- subject level --- 
figure
set(gcf,'Position',[100 100 800 200])

% --- Plot precue T1 ---
subplot 151
val = circ_mean(fittedP.linear2Hz.valR1_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(1,:),'EdgeColor',p.cueColors(1,:)); % fitted phase histogram
hold on 
[mu ul ll] = circ_mean(val); 
r = circ_r(val); 
polarplot([mu mu],[0 r],'Color','b','LineWidth',2)
radius = 3;
th = ll:0.01:ul;
polarplot(th,radius+zeros(size(th)),'Color',[0.5 0.5 0.5],'LineWidth',2)
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T1')

% --- Plot precue T2 ---
subplot 152
val = circ_mean(fittedP.linear2Hz.valR2_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(2,:),'EdgeColor',p.cueColors(2,:)); % fitted phase histogram
hold on 
[mu ul ll] = circ_mean(val); 
r = circ_r(val); 
polarplot([mu mu],[0 r],'Color','r','LineWidth',2) 
th = ll:0.01:ul;
polarplot(th,radius+zeros(size(th)),'Color',[0.5 0.5 0.5],'LineWidth',2)
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T2')

% --- Plot precue T1 - precue T2 difference -> averaged  ---
subplot 153
val = circ_mean(fittedP.linear2Hz.valR_diffMin,[],2); % average sessions  
polarhistogram(val,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1 - T2\n avg diff'))

% --- Plot precue T1, precue T2 averaged -> difference --- 
subplot 154
val1 = circ_mean(fittedP.linear2Hz.valR1_min,[],2);
val2 = circ_mean(fittedP.linear2Hz.valR2_min,[],2);
val = circ_dist(val1,val2);
fittedP.linear2Hz.valR_avgThenDiff = val; 
polarhistogram(val,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1, T2\n avg, then diff'))

% --- Plot precue T1 - precue T2 abs difference ---
subplot 155
val = circ_mean(fittedP.linear2Hz.valR_diffMinAbs,[],2); % average sessions  
polarhistogram(val,edges); % fitted phase histogram
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('abs(Precue T1 - T2)')

sgtitle(sprintf('Subjects (n = %d)',...
    nSubs/2))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_PermutationsHistogram_sessions%d',...
        mdlFitType,...
        nSubs/2);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Plot precue T1 and precue T2 phase distribution 95% CI same fig 
figure
set(gcf,'Position',[100 100 230 170])

% --- Plot precue T1 ---
val = circ_mean(fittedP.linear2Hz.valR1_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(1,:),'EdgeColor',p.cueColors(1,:),'FaceAlpha',1); % fitted phase histogram
hold on 
[mu ul ll] = circ_mean(val); 
r = circ_r(val); 
polarplot([mu mu],[0 r],'Color',p.darkColors(1,:),'LineWidth',2)
radius = 3;
th = ll:0.01:ul;
polarplot(th,radius+zeros(size(th)),'Color',p.darkColors(1,:),'LineWidth',2)
pax = gca;
pax.ThetaAxisUnits = 'radians';

% --- Plot precue T2 ---
val = circ_mean(fittedP.linear2Hz.valR2_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(2,:),'EdgeColor',p.cueColors(2,:),'FaceAlpha',1); % fitted phase histogram
hold on 
[mu ul ll] = circ_mean(val); 
r = circ_r(val); 
polarplot([mu mu],[0 r],'Color',p.darkColors(2,:),'LineWidth',2) 
th = ll:0.01:ul;
polarplot(th,radius+zeros(size(th)),'Color',p.darkColors(2,:),'LineWidth',2)
pax = gca;
pax.ThetaAxisUnits = 'radians';

pax.FontSize = 16;
pax.FontName = 'Helvetica-Light';
rticks([1,2,3])
thetaticks(0:pi/4:2*pi)

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_Summary_sessions%d_PrecueT1T2',...
        mdlFitType,...
        nSubs/2);
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Plot precue T1 and precue T2 DIFFERENCE phase distribution 95% CI same fig 
figure
set(gcf,'Position',[100 100 230 170])

% --- Plot precue T1 - precue T2 ---
val = fittedP.linear2Hz.valR_avgThenDiff; % average sessions  
polarhistogram(val,edges,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7],'FaceAlpha',0.7); % fitted phase histogram
hold on 
[mu ul ll] = circ_mean(val); 
r = circ_r(val); 
polarplot([mu mu],[0 r],'Color','k','LineWidth',2)
radius = 2;
th = ll:0.01:ul;
polarplot(th,radius+zeros(size(th)),'Color','k','LineWidth',2)
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.FontSize = 16;
pax.FontName = 'Helvetica-Light';
rticks([1,2])
thetaticks(0:pi/4:2*pi)

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_Summary_sessions%d_PrecueT1T2Diff',...
        mdlFitType,...
        nSubs/2);
    saveas(gcf,sprintf('%s/%s.svg', figDir, figTitle))
end

%% Plot group (phase) by amplitude --- Subjects --- 
figure
set(gcf,'Position',[100 100 800 200])

% --- Plot precue T1 ---
subplot 161
% rho is amplitude 
rho = mean(fittedP.linear2Hz.valAmp1_min,2); 
% theta is phase 
theta = circ_mean(fittedP.linear2Hz.valR1_min,[],2);
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T1')

% --- Plot precue T2 ---
subplot 162
% rho is amplitude 
rho = mean(fittedP.linear2Hz.valAmp2_min,2); 
% theta is phase 
theta = circ_mean(fittedP.linear2Hz.valR2_min,[],2);
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title('Precue T2')

% --- Plot precue T1 - precue T2 difference ---
subplot 163
% rho is amplitude difference of precue T1 and T2 
rho = mean(fittedP.linear2Hz.valAmp1_min,2) - mean(fittedP.linear2Hz.valAmp2_min,2); 
% theta is phase 
theta = circ_mean(fittedP.linear2Hz.valR_diffMin,[],2); 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Diff->Avg: Precue T1 - T2\nAmp difference'))

% --- Plot precue T1 - precue T2 difference ---
subplot 164
% rho is amplitude mean of precue T1 and T2 
rho = (mean(fittedP.linear2Hz.valAmp1_min,2)+mean(fittedP.linear2Hz.valAmp2_min,2))/2; 
% theta is phase 
theta = circ_mean(fittedP.linear2Hz.valR_diffMin,[],2); 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Diff->Avg: Precue T1 - T2\nAmp mean'))

% --- Plot precue T1 - precue T2 difference ---
subplot 165
% rho is amplitude difference of precue T1 and T2 
rho = mean(fittedP.linear2Hz.valAmp1_min,2) - mean(fittedP.linear2Hz.valAmp2_min,2); 
% theta is phase 
theta = fittedP.linear2Hz.valR_avgThenDiff; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Avg->Diff: Precue T1 - T2\nAmp difference'))

% --- Plot precue T1 - precue T2 difference ---
subplot 166
% rho is amplitude mean of precue T1 and T2 
rho = (mean(fittedP.linear2Hz.valAmp1_min,2)+mean(fittedP.linear2Hz.valAmp2_min,2))/2; 
% theta is phase 
theta = fittedP.linear2Hz.valR_avgThenDiff; 
polarscatter(theta,rho,'Marker','+'); 
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Avg->Diff: Precue T1 - T2\nAmp mean'))

sgtitle(sprintf('Group (n = %d subjects) Phase x Amplitude',...
    nSubs/2))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear%dHz_Phase-Amp_Subjects%d',...
        mdlFitType,freq,...
        nSubs/2);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Plot group (phase) --- subject level (OVERLAID) --- 
figure
set(gcf,'Position',[100 100 300 200])

val = circ_mean(fittedP.linear2Hz.valR1_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(1,:),'FaceAlpha',1); % fitted phase histogram
hold on
val = circ_mean(fittedP.linear2Hz.valR2_min,[],2); % average sessions  
polarhistogram(val,edges,'FaceColor',p.cueColors(2,:),'FaceAlpha',0.7); 
pax = gca;
pax.ThetaAxisUnits = 'radians';

sgtitle(sprintf('Group averaged to subjects (n = 10)'))

if saveFigs
    figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_%s_Linear2Hz_Phase_PermutationsHistogram_subjects_overlaid',...
        mdlFitType);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Phase stats on session and subjects
clear stats
levels = {'sessions','subjects'};
% idx = [find(contains(fields,'valR1_min'));find(contains(fields,'valR2_min'));...
%     find(contains(fields,'valR_diffMin'))]'; 
for iL = 1:numel(levels)
    if strcmp(levels{iL},'sessions')
        level = 'sessions';
        fields = fieldnames(group); 
        idx = [find(contains(fields,'valR1_min'));find(contains(fields,'valR2_min'));...
            find(contains(fields,'valR_diffMin'))]';
    elseif strcmp(levels{iL},'subjects')
        level = 'subjects';
        fields = fieldnames(fittedP.linear2Hz); 
        idx = [find(contains(fields,'valR1_min'));find(contains(fields,'valR2_min'));...
            find(contains(fields,'valR_avgThenDiff'));find(contains(fields,'valR_diffMin'))]';
    else
        error('level for statistical test not found')
    end

    for iV = idx
        clear val
        if strcmp(levels{iL},'sessions')
            level = 'sessions';
            val = group.(fields{iV})';
        elseif strcmp(levels{iL},'subjects')
            level = 'subjects';
            if strcmp(fields{iV},'valR_avgThenDiff')
                val = fittedP.linear2Hz.(fields{iV}); 
            else
                val = circ_mean(fittedP.linear2Hz.(fields{iV}),[],2);
            end
        else
            error('level for statistical test not found')
        end
        % --- Mean + CI ---
        [mu,ul,ll] = circ_mean(val,[],1);
        stats.(level).(fields{iV}).mu = mu;
        stats.(level).(fields{iV}).ul = ul;
        stats.(level).(fields{iV}).ll = ll;

        % --- Decriptive stats --- 
        descriptive = circ_stats(val);
        kappa = circ_kappa(val);
        stats.(level).(fields{iV}).stats = descriptive; 
        stats.(level).(fields{iV}).stats.kappa = kappa; 

        % --- Rayleighs z test ---
        testType = 'rayleighsZ'; % tests for uniformity   
        [pval,z] = circ_rtest(val);
        stats.(level).(fields{iV}).(testType).circ_test_name = 'circ_rtest'; 
        stats.(level).(fields{iV}).(testType).level = level; 
        stats.(level).(fields{iV}).(testType).pval = pval;
        stats.(level).(fields{iV}).(testType).z = z;
        
        % --- Parametric Watson-Williams multi-sample test --- 
        % H0: the populations have equal means
        % HA: the populations have unequal means 
        clear val_precueT1 val_precueT2 
        val_precueT1 = group.valR1_min';
        val_precueT2 = group.valR2_min';
        testType = 'watsonwilliams';
        [pval, table] = circ_wwtest(val_precueT1,val_precueT2);
        stats.(level).(fields{iV}).(testType).circ_test_name = 'circ_wwtest'; 
        stats.(level).(fields{iV}).(testType).level = level; 
        stats.(level).(fields{iV}).(testType).pval = pval;
        stats.(level).(fields{iV}).(testType).table = table; 

       % --- Parametric 2-way ANOVA for circular data with interaction --- 
       % Independent samples t-test
       testType = 'twowayANOVA'; 
       clear val 
       val = [val_precueT1;val_precueT2]; 
       n = size(val); 
       idp = ones(n); % factor 1: precue 
       idp(1:(n/2)) = 1; 
       idp((n/2)+1:n) = 2; 
       idq = ones(n); % factor 2: session 
       idq(2:2:40) = 2; 
       inter = 1; % logical, report interaction 
       fn = {'precue','session'};
       [pval table] = circ_hktest(val, idp, idq, inter, fn);
       stats.(level).(fields{iV}).(testType).circ_test_name = 'circ_hktest';
       stats.(level).(fields{iV}).(testType).level = level;
       stats.(level).(fields{iV}).(testType).pval = pval;
       stats.(level).(fields{iV}).(testType).table = table;
        
    end
end
group.stats = stats; 

%% Prepare for export 
% Long format data for R
exportCSV = 0; 
if exportCSV
    fitTypes = {'linear','linear2Hz'};
    for iF = 2 % just 'linear2Hz' 1:numel(fitTypes)
        clear V
        fields = fieldnames(fittedP.(fitTypes{iF}));
        tableHeaders = {'subject','session',...
            'valR1_min','valT1_min',...
            'valR2_min','valT2_min',...
            'valR_diffMin','valT_diffMin',...
            'valR_diffMinAbs','valT_diffMinAbs'};
        count = 1;
        for iS = 1:10 % subjects
            for iSession = 1:2
                for iV = 1:numel(fields)
                    V.subject(count) = iS;
                    V.session(count) = iSession;
                    V.(tableHeaders{iV+2})(count) = fittedP.(fitTypes{iF}).(fields{iV})(iS,iSession);
                end
                count = count+1;
            end
        end

        % --- Export to csv for R ANOVA ---
        sz = [numel(V.(tableHeaders{1})) numel(tableHeaders)];

        varTypes = repmat("double",[1 numel(tableHeaders)]);
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

        for iT = 1:numel(tableHeaders)
            T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
        end

        dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
        csvName = sprintf('TANoise_ITPCFit_%s_Phase_%s_%s',fitTypes{iF},mdlFitType,dateStr);
        csvPath = sprintf('%s/%s.csv', csvDir, csvName);
        writetable(T,csvPath)
    end
end





