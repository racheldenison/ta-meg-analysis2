%% ITPC temporal expectation fit by precue 
% fmincon requires the optimization toolbox 
% set path to high level ta-meg-analysis2
% Fit ITPC data, separately to precue T1 and T2 
% August 2023 Karen 

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian (persona), karen (lab)

addpath(genpath(pwd))

freq = 3; 
dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPCfit_separate1_%dHz/%s',user,freq,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Prepare data 
% --- Load ITPC data ---
load(sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/unused/groupA_ITPCspectrogram_byAtt.mat',user)) % A variable 
% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

% --- Data settings ---
% sampling = 1:1:7001; % 10 ms 
foi = 20; % frequency of interest, Hz
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 

Fs = 1000; % sampling frequency 

fitLevel = 'session'; 
baselineToi = abs(p.tstart)+p.eventTimes(1)-100:abs(p.tstart)+p.eventTimes(1); 

% -- MEG settings --- 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

% --- Extract data into precue conds at desired foi --- 
cueNames = {'all','cueT1','cueT2'};
fitLevels = {'session','subject','group'}; 
clear data 
for i = 1:numel(cueNames)
    data.(cueNames{i}).session = squeeze(A.(cueNames{i}).session(foi,:,:)); 
    data.(cueNames{i}).subject = squeeze(A.(cueNames{i}).subject(foi,:,:)); % frequency x time x subjects
    data.(cueNames{i}).group   = squeeze(mean(data.(cueNames{i}).subject,2)); % average across subjects
end

%% --- Fit settings --- 
paramNames = {'intercept','slope','amplitude','phase'};
nVars = numel(paramNames); 

% === SETTINGS ===
fitLevel = 'session'; % 'session' 'group'
cueLevel = {'cueT1','cueT2'};
fitTypes = {'linear2Hz'}; % 'linear2Hz' 'linear'
o_bads = 1; % if 1, uses bads, otherwise uses fmincon
randCoeffs = 1; % uses random starting coefficients from the search space, otherwise uses hardcoded intial coefs
% ================

clear mdlFit lb ub 
mdlFit.A = []; % equality contraints
mdlFit.b = []; 
mdlFit.Aeq = []; % inequality contraints
mdlFit.beq = []; 

% === Define upper and lower bounds ===

% --- Intercept ---
idx = find(contains(paramNames,'intercept')); 
lb(idx) = 0; 
ub(idx) = 1; 

% --- Slope --- 
idx = find(contains(paramNames,'slope')); 
lb(idx) = -2; 
ub(idx) = 2; % 0.001; % arbitrarily large, for the scaling? 

% --- Amplitude --- 
idx = find(contains(paramNames,'amplitude')); 
lb(idx) = 0; 
ub(idx) = 1; 

% --- Phase shift ---
idx = find(contains(paramNames,'phase')); 
lb(idx) = 0; 
ub(idx) = 2*pi; % 500 
if o_bads
    optionsBads = bads('defaults');
    optionsBads.PeriodicVars = idx;
    lb(idx) = -pi; 
    ub(idx) = pi;
end

% --- Frequency --- (fix to 2 Hz?) 
% idx = find(contains(paramNames,'freq')); 
% lb(idx) = 2; % *pi/(500); % 1.5 Hz 1.7
% ub(idx) = 2; % *pi/(500); % 0.0016; % 2.5 Hz  2.3

mdlFit.lb = lb; 
mdlFit.ub = ub; 

mdlFit.hardlb = mdlFit.lb;
mdlFit.hardub = mdlFit.ub;

% --- Nonlinear constraints --- 
nonlcon = []; 

% --- Options --- 
options = optimoptions('fmincon','Display','iter'); 

% if o_bads
%     optionsBads = bads('defaults');
%     optionsBads.PeriodicVars = 4;
%     mdlFit.lb(4) = -2.5;
%     mdlFit.ub(4) = 2.5;
%     mdlFit.hardlb = mdlFit.lb;
%     mdlFit.hardub = mdlFit.ub;
% end

%% --- Define initial search array --- 
nGrain = 100; 
clear x0s
x0s = NaN([nVars,nGrain]);
paramNames = {'intercept','slope','amplitude','phase'};
for i = 1:nVars
    % if strcmp(paramNames{i},'slope')
    %     x0s(i,:) = linspace(lb(i),0,nGrain); % slope max to 0? 
    % else
        x0s(i,:) = linspace(mdlFit.lb(i),mdlFit.ub(i),nGrain); % initial coefficients for search
    % end
end

%% --- Fit model ---
clear x0 x0Perm
nPerms = 100; % 100
for iS = 1:size(data.all.(fitLevel),2) % subjects
    for iP = 1:nPerms
        txt = sprintf('Fitting %d of %d, permutation %d of %d...',iS,size(data.all.(fitLevel),2),iP,nPerms); 
        disp(txt)
        % --- Randomly pick starting coefficients from search grid ---
        if randCoeffs
            for iV = 1:nVars
                vals = randperm(nGrain);
                idx(iV) = vals(1);
                x0(iV) = x0s(iV,idx(iV));
            end
        else
            % --- Hard code starting coefs ---- EDIT ----
            x0 = [0.3 2 0.1 0];
        end

        for iC = 1:numel(cueLevel)
            clear dataFit
            dataFit = data.(cueLevel{iC}).(fitLevel)(tIdx,iS)';

            for iF = 1:numel(fitTypes)
                % --- Save fit settings ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).x0(iP,iS,:) = x0; % save starting coefficients
                if o_bads
                    algo = 'bads';
                else
                    algo = 'fmincon';
                end
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).algo = algo; % save optimization algorithm

                switch fitTypes{iF}
                    case 'linear2Hz'
                        disp('Fitting linear + 2Hz ...')
                        paramNames = {'intercept','slope','amplitude','phase'};
                        % --- Define objective function ---
                        clear fun
                        fun = @(x)meg_objectiveFunction1(x,dataFit,t,Fs,paramNames,fitTypes{iF},freq);
                        % --- Do fit ---
                        if o_bads
                            [solution,fval,exitflag,output] = bads(fun, x0, mdlFit.hardlb, mdlFit.hardub, mdlFit.lb, mdlFit.ub,[],optionsBads);
                        else
                            [solution,fval,exitflag,output] = fmincon(fun, x0, mdlFit.A, mdlFit.b, mdlFit.Aeq, mdlFit.beq, mdlFit.lb, mdlFit.ub);
                        end
                    case 'linear'
                        disp('Fitting linear ...')
                        x0L = x0(1:2);
                        paramNames = {'intercept','slope'};
                        % --- Define objective function ---
                        clear fun
                        fun = @(x)meg_objectiveFunction1(x,dataFit,t,Fs,paramNames,fitTypes{iF});
                        idx = 1:numel(paramNames);
                        % --- Do fit ---
                        if o_bads
                            [solution,fval,exitflag,output] = bads(fun, x0L, [], [], mdlFit.lb(1:2), mdlFit.ub(1:2));
                        else
                            [solution,fval,exitflag,output] = fmincon(fun, x0L, mdlFit.A, mdlFit.b, mdlFit.Aeq, mdlFit.beq, mdlFit.lb(idx), mdlFit.ub(idx));
                        end
                end
                % --- Save parameter names ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames = paramNames; 
                % -- Save fixed freq --- 
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fixedFreq = freq; % hz 
                % --- Save fitted parameters ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iS,:) = solution;
                % --- Save error ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(iP,iS) = fval;
                % --- Save predicted y ---
                % clear E yhat
                % [E,yhat] = meg_objectiveFunction1(solution,dataFit,t,Fs,paramNames,fitTypes{iF});
                % mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(iP,iS,:) = yhat;
            end
        end
    end
    % Save workspace and model fit
    filename = sprintf('ModelFit_SeparatePrecueT1T2_%dHz_%s',freq,dateStr);
    save(filename,'-v7.3')
end

%% Save additional info? -- should be stored now, temp solve --- 
% for iS = 1:size(data.all.(fitLevel),2) % subjects
%     for iP = 1:nPerms
%         for iC = 1:numel(cueLevel)
%             for iF = 1:numel(fitTypes)
%                 switch fitTypes{iF}
%                     case 'linear2Hz'
%                         paramNames = {'intercept','slope','amplitude','phase'};
%                     case 'linear'
%                         x0L = x0(1:2);
%                         paramNames = {'intercept','slope'};
%                 end
%                 % --- Save parameter names ---
%                 mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames = paramNames; 
%                 % -- Save fixed freq --- 
%                 mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fixedFreq = freq; % hz 
%             end
%         end
%     end
% end
% % Save workspace and model fit
% filename = sprintf('ModelFit_SeparatePrecueT1T2_%dHz_%s',freq,dateStr);
% save(filename,'-v7.3')

%% --- Figure (data + fit) ---
for iF = 1:numel(fitTypes)
    for iS = 1:size(data.all.(fitLevel),2) % sessions
    figure
    switch fitLevel
        case 'session'
            set(gcf,'Position',[100 100 450 300])
        case 'subject'
            set(gcf,'Position',[100 100 450 300])
        otherwise
            error('Specify session- or subject-level fit')
    end
    hold on
    
    for iC = 1:numel(cueLevel)
        % --- Data (by precue) ---
        dataFit = data.(cueLevel{iC}).(fitLevel)(:,iS)';

        % --- Plot data ---
        plot(p.t,dataFit,'LineWidth',1,'Color',p.cueColors(iC,:))

        % --- Plot model fit ---
        [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS)  );
        if iC == 1
            idx1 = idx; 
        elseif iC ==2 
            idx2 = idx; 
        end
        fittedX = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:);
        switch fitTypes{iF}
            case 'linear2Hz'
                paramNames = {'intercept','slope','amplitude','phase'};
            case 'linear'
                paramNames = {'intercept','slope'};
        end
        [~,yhat] = meg_objectiveFunction1(fittedX,dataFit(toi),t,Fs,paramNames,fitTypes{iF},freq);
        fitColors = {'b','r'};
        plot(t,yhat,'--','LineWidth',1,'Color',fitColors{iC})
    end

    % --- Plot event lines ---
    for i = 1:numel(p.eventTimes)
        xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
    end

    % --- Format ---
    meg_figureStyle
    xlim([-100 2400])
    xlabel('Time (ms)')
    ylabel('ITPC')

    % --- Titles --
    switch fitTypes{iF}
        case 'linear2Hz'
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\\pi\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\\pi\nfreq = %0.2f',...
                mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2),...
                mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,3),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,4)/pi,...
                mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,2),...
                mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,3),mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,4)/pi,...
                freq);
        case 'linear'
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f\nintercept2 = %0.2f, slope2 = %0.3f',...
                mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2),...
                mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS,2));
    end
    titleText1 = sprintf('Separate precue T1 and T2 fits (%s), fval = %0.2e, %s',...
        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
        mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(iS),...
        und2space(sessionNames{iS})); 
    titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted); 
    title(titleText)
    
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.6;
    ax.TitleFontWeight = 'normal';

    % --- Save fig ---
    dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
    figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_%s_%s',sessionNames{iS},fitTypes{iF},dateStr);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    end

    % lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast');
    % lgd.FontSize = 10;
end

%% Polar histogram (phase summary subject permutations) ********* 
% load('ModelFit_SeparatePrecueT1T2_2Hz_231006.mat') 

% load('ModelFit_SeparatePrecueT1T2_1Hz_231014.mat')
% load('ModelFit_SeparatePrecueT1T2_3Hz_231014.mat')
% load('ModelFit_SeparatePrecueT1T2_4Hz_231014.mat')
% load('ModelFit_SeparatePrecueT1T2_5Hz_231014.mat')
% load('ModelFit_SeparatePrecueT1T2_6Hz_231014.mat')

mdlFitType = 'separate'; 
[group,fittedP,stats] = meg_plotPhaseHistogram(mdlFit,paramNames,figDir,mdlFitType); 

%% Save min solution 
for iF = 1:numel(fitTypes)
    for iS = 1:size(data.all.(fitLevel),2) % sessions 
        for iC = 1:numel(cueLevel)
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS));
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minFval(iS) = minVal;
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:) = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:); 
        end
    end
end

% Average min fval to subjects 
for iF = 1:numel(fitTypes)
    for iC = 1:numel(cueLevel)
    end
end

% Average to subjects 
for iF = 1:numel(fitTypes)
    for iC = 1:numel(cueLevel)
        subCount = 0;
        for iS = 1:2:size(data.all.(fitLevel),2) % sessions
            subCount = subCount+1;
            clear vals1 vals2 vals12 meanVals12
            vals1 = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:);
            vals2 = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS+1,:);
            vals12 = cat(1,vals1,vals2);
            meanVals12 = mean(vals12,1,'omitnan');
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolutionSub(subCount,:) = meanVals12;
        end
    end
end

%% --- Plot group params by precue --- 231013
% fig styling 
sz = 100; 
for iF = 2 % 1:numel(fitTypes)
    switch fitTypes{iF}
        case 'linear2Hz'
            disp('Plotting linear + 2Hz ...')
            paramNames = {'intercept','slope','amplitude','phase'};
    end
    % --- group --- 
    figure
    set(gcf,'Position',[100 100 300 500])
    for iV = 1:numel(paramNames)
        for iC = 1:numel(cueLevel)
            val = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolutionSub(:,iV); 
            meanVal = mean(val,1); % average subjects
            err = std(val)/sqrt(numel(val)); 
            subplot (numel(paramNames),1,iV)
            hold on
            figureStyle
            scatter(iC,meanVal,sz,'MarkerFaceColor',p.cueColors(iC,:))
            er = errorbar(iC,meanVal,err);
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            er.LineWidth = 2; 
            er.CapSize = 0; 
            xticks([1 2])
            xticklabels({'Precue T1','Precue T2'})
            xlim([0 3])
            ylabel(sprintf(paramNames{iV}))
        end
    end
end

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_Separate_%s_%s_Avg_ByPrecue_n10',fitTypes{iF},dateStr);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% --- Plot precue T1 x precue T2 unity plots 231013 --- 
for iF = 2 % 1:numel(fitTypes)
    switch fitTypes{iF}
        case 'linear2Hz'
            disp('Plotting linear + 2Hz ...')
            paramNames = {'intercept','slope','amplitude','phase'};
    end
    % --- group ---
    figure
    set(gcf,'Position',[100 100 300 500])
    for iV = 1:numel(paramNames)
        % --- Precue T1 --- 
        val1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).minSolutionSub(:,iV);
        meanVal1 = mean(val1,1); % average subjects
        % --- Precue T2 --- 
        val2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).minSolutionSub(:,iV);
        meanVal2 = mean(val1,2); % average subjects
        subplot (numel(paramNames),1,iV)

        hold on
        figureStyle
        scatter(val1,val2,sz,'MarkerFaceColor',p.cueColors(iC,:))
        refline(1,0)
        title(sprintf(paramNames{iV}))
        xlabel('Precue T1')
        ylabel('Precue T2')
        axis square
    end
end

%% --- Correlate amplitude and phase (session) --- 
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
polarhistogram(group.valR2_min,edges); % fitted phas histogram
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

%% Polar histogram (phase summary subject permutations) 
edges = 0:pi/24:2*pi; 
for iS = 1:size(mdlFit.linear2Hz.cueT1.session.fval,2) % :size(data.all.(fitLevel),2) % subjects
    figure
    subplot 131
    % --- Calculate phase difference ---
    % valT_diff = abs(mdlFit.linear2Hz.(cueLevel{1}).(fitLevel).solution(idx,iS,4) - mdlFit.linear2Hz.(cueLevel{2}).(fitLevel).solution(idx,iS,4));
    % valR_diff = t2rad(valT_diff);
    valR_diff = circ_dist(t2rad(mdlFit.linear2Hz.(cueLevel{1}).(fitLevel).solution(:,iS,4)), t2rad(mdlFit.linear2Hz.(cueLevel{2}).(fitLevel).solution(:,iS,4)));
    % --- Plot fitted phases (radians) Precue T1 ---
    for iC = 1
        rho = 1;
        valT = mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution(:,iS,4);
        valRad = t2rad(valT);
        polarhistogram(valRad,edges)
        hold on
        txt = sprintf('%0.2f deg = %0.2f rad',valT,valRad);
    end
    % --- Plot fitted phase w min fval - Precue T1 ---
    [minVal,idx1] = min(  mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(:,iS)  );
    valT = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx1,iS,4);
    valRad = t2rad(valT);
    polarplot([0 valRad],[0 100],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(sprintf('Fitted phase (rad)\nPrecue T1'))

    subplot 132
    % --- Plot fitted phases (radians) Precue T2 ---
    for iC = 2
        rho = 1;
        valT = mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution(:,iS,4);
        valRad = t2rad(valT);
        polarhistogram(valRad,edges)
        hold on 
        txt = sprintf('%0.2f deg = %0.2f rad',valT,valRad);
    end
    % --- Plot fitted phase w min fval - Precue T2 ---
    [minVal,idx2] = min(  mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).fval(:,iS)  );
    valT = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx2,iS,4);
    valRad = t2rad(valT);
    polarplot([0 valRad],[0 100],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(sprintf('Fitted phase (rad)\nPrecue T2'))

    subplot 133
    % --- Plot phase difference ---
    polarhistogram(valR_diff,edges)
    hold on
    % --- Plot fitted phase w min fval - Precue T1-Precue T2 ---
    valR_diff = circ_dist(t2rad(mdlFit.linear2Hz.(cueLevel{1}).(fitLevel).solution(idx1,iS,4)), t2rad(mdlFit.linear2Hz.(cueLevel{2}).(fitLevel).solution(idx2,iS,4)));
    polarplot([0 valR_diff],[0 100],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(sprintf('Fitted phase (rad) difference:\nPrecue T1 - Precue T2'))

    sgtitle(sprintf('%s, nPerm = %d',...
        und2space(sessionNames{iS}),...
        size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,1)))

    % --- Save fig ---
    figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Phase_PermutationsHistogram_%d',...
        sessionNames{iS},...
        size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,1));
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Polar histogram (phase summary group permutations) 
edges = 0:pi/24:2*pi; 
freq = 2; % Hz 
% --- Get phase with min fval ---
clear minVal1 idx1 minVal2 idx2 valT_precueT1 valT_precueT2 valR_precueT1 valR_precueT2 valR_diff
for iS = 1:size(mdlFit.linear2Hz.cueT1.session.fval,2)
    [minVal1(iS),idx1(iS)] = min(  mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(:,iS)  );
    [minVal2(iS),idx2(iS)] = min(  mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).fval(:,iS)  );

    % --- Phase in t ---
    valT_precueT1(iS,:) = mdlFit.linear2Hz.(cueLevel{1}).(fitLevel).solution(:,iS,find(contains(paramNames,'phase')));
    valT_precueT2(iS,:) = mdlFit.linear2Hz.(cueLevel{2}).(fitLevel).solution(:,iS,find(contains(paramNames,'phase')));
end
% --- Phase in rad ---
valR_precueT1 = t2rad( valT_precueT1, freq );
valR_precueT2 = t2rad( valT_precueT2, freq );
% --- Calculate phase difference ---
valR_diff = circ_dist(valR_precueT1,valR_precueT2);

% -- Get min phase values --- 
for iS = 1:size(mdlFit.linear2Hz.cueT1.session.fval,2)
    valT_precueT1_min(iS) = valT_precueT1(iS,idx1(iS));
    valT_precueT2_min(iS) = valT_precueT2(iS,idx2(iS));
end
% --- Phase in rad ---
valR_precueT1_min = t2rad( valT_precueT1_min, freq );
valR_precueT2_min = t2rad( valT_precueT2_min, freq );
% --- Calculate phase difference ---
valR_diff_min = circ_dist(valR_precueT1_min,valR_precueT2_min); 

figure
% --- Plot group polar histogram ---
subplot 131
% --- Plot fitted phases (radians) Precue T1 ---
rho = 1;
polarhistogram(valR_precueT1,edges)
hold on
% --- Plot fitted phase w min fval - Precue T1 ---
for iS = 1:size(mdlFit.linear2Hz.cueT1.session.fval,2)
    polarplot([1 valR_precueT1_min(iS)],[0 300],'r')
end
% --- Formatting ---
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1'))

subplot 132
% --- Plot fitted phases (radians) Precue T2 ---
rho = 1;
polarhistogram(valR_precueT2,edges)
hold on
% --- Plot fitted phase w min fval - Precue T2 ---
for iS = 1:size(mdlFit.linear2Hz.cueT2.session.fval,2)
    polarplot([1 valR_precueT2_min(iS)],[0 300],'r')
end
% --- Formatting ---
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T2'))

subplot 133
% --- Plot fitted phases (radians) Difference ---
rho = 1;
polarhistogram(valR_diff,edges)
hold on
% --- Plot fitted phase w min fval - Precue T2 ---
for iS = 1:size(mdlFit.linear2Hz.cueT2.session.fval,2)
    polarplot([1 valR_diff_min(iS)],[0 300],'r')
end
% --- Formatting ---
pax = gca;
pax.ThetaAxisUnits = 'radians';
title(sprintf('Precue T1 - Precue T2'))

sgtitle(sprintf('Fitted phase (rad), n = %d, nPerm = %d',...
    size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,2),...
    size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,1)))

% --- Save fig ---
figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Phase_PermutationsHistogram_nPerm%d_n%d',...
    size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,1),...
    size(mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).solution,2));
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% Session specific model comparison? 
for iS = 1:2 % size(mdlFit.linear2Hz.cueT1.(fitLevel).fval,2)
    for iC = 1:numel(cueLevel)
        % --- Find parameters resulting in min fval (linear only) ---
        clear minVal idx
        [minVal,idx] = min(  mdlFit.linear.(cueLevel{iC}).(fitLevel).fval(:,iS)  );
        sse1 = mdlFit.linear.(cueLevel{iC}).(fitLevel).fval(idx,iS);
        k1 = 2; 
        
        % --- Find parameters resulting in min fval (linear + 2Hz) ---
        clear minVal idx
        [minVal,idx] = min(  mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).fval(:,iS)  );
        sse2 = mdlFit.linear2Hz.(cueLevel{iC}).(fitLevel).fval(idx,iS);
        k2 = 4; 
        
        % --- Model comparison ---
        n = 20; % How to determine? 
        fval(iS,iC) = meg_nestedModelComparison(sse1,sse2,k1,k2,n);
    end
end

%% Get only min fval solution (current) 
for iF = 1:numel(fitTypes)
    for iS = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution,2)
        for iC = 1:2
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:,iS)  );
            fittedX = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,iS,:);
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:) = fittedX; 
        end
    end
end

%% Prepare data for figures of fitted params
% Restructure data into 10 sub x 2 sessions x 2 precues
clear fittedP
for iF = 1:numel(fitTypes)
    for iV = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).minSolution,2) % params
        subCount = 0;
        for iS = 1:2:20 % sessions
            subCount = subCount+1;
            clear val valCT1 valCT2
            valCT1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).minSolution(iS:iS+1,iV);
            valCT2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).minSolution(iS:iS+1,iV);
            val = squeeze(cat(4,valCT1,valCT2));
            fittedP.(fitTypes{iF}).(paramNames{iV})(subCount,:,:) = val;
        end
        % --- Transform phase from t to rad space ---
        % if strcmp(paramNames{iV},'phase')
        %     clear val 
        %     val = fittedP.(fitTypes{iF}).(paramNames{iV}); 
        %     % val = (val*2*pi)/Fs; % if -500 to 500 bounds
        %     val = (val*4*pi)/Fs;
        %     fittedP.(fitTypes{iF}).phaserad = val; 
        % end
    end
end

% --- Compile data and fit baseline ---
for iF = 1:numel(fitTypes)
    subCount = 0;
    for iS = 1:2:20 % sessions
        subCount = subCount+1;
         
        % --- Data ---
        clear val valCT1 valCT2
        valCT1 = mean(data.(cueLevel{1}).(fitLevel)(baselineToi,iS:iS+1),1,'omitnan');
        valCT2 = mean(data.(cueLevel{2}).(fitLevel)(baselineToi,iS:iS+1),1,'omitnan');

        val = squeeze(cat(4,valCT1,valCT2));
        fittedP.(fitTypes{iF}).baseline_data(subCount,:,:) = val;

        % --- Fit ---
        clear val valCT1 valCT2
        valCT1 = mean(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).yhat(iS:iS+1,baselineToi),2,'omitnan');
        valCT2 = mean(mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).yhat(iS:iS+1,baselineToi),2,'omitnan');

        val = squeeze(cat(4,valCT1,valCT2));
        fittedP.(fitTypes{iF}).baseline_yhat(subCount,:,:) = val;
    end
end

%% --- Plot fitted params ---
for iF = 1:numel(fitTypes)
    meg_plotFittedParam1(  fittedP.(fitTypes{iF}), p, fitTypes{iF}, figDir) % 10 sub x 2 sessions x 2 precues
end

%% Ttest on phase rad difference between precue T1 and precue T2 
% see meg_phaseDiff 

%% --- ttests on fitted model params (Precue T1 vs Precue T2) --- 
% paramNamesAdd = paramNames; 
% paramNamesAdd{end+1} = 'phaserad';
% paramNamesAdd{end+1} = 't0'; 
% 
% clear fittedVal
% for iF = 1:numel(fitTypes)
%     for iV = 1:numel(paramNamesAdd)
%         for iC = 1:numel(cueLevel)
%             for iS = 1:size(data.all.(fitLevel),2)
%                 clear val
%                 if iV<=5
%                     % Average to sessions ?
%                     val = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS);
%                     fittedVal(iV,iC,iS,iF) = val;
%                 end
%             end
%         end
% 
%         % Convert phase from t to rad, append
%         switch paramNamesAdd{iV}
%             case {'phase'} % circular stats
%                 fittedVal = deg2rad(fittedVal);
%                 val = (val*2*pi)/Fs; % covert from t to rad
%                 % fittedVal(iV,iC,iS,iF) = (fittedVal*2*pi)/Fs;
%                 % fittedVal(iV,iC,iS,iF) = val; 
%                 mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS) = fittedVal(iV,iC,iS,iF);
%             % case {'intercept'}
%             %     fittedX = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(1,:,iS);
%             %     dataFit = data.(cueLevel{iC}).(fitLevel);
%             %     dataFit = dataFit(:,iS);
%             %     [E,yhat] = objectiveFunction(fittedX,dataFit,t,Fs);
%             %     val = yhat(1); 
%             %     fittedVal(iV,iC,iS,iF) = val; 
%             %     mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS) = fittedVal(iV,iC,iS,iF);
%         end
%     end
% end

%% Compare fit errors of linear vs linear+2Hz
figure 
for iC = 1:numel(cueLevel)
    subplot (2,1,iC)
    hold on 
    meg_figureStyle
    for iF = 1:numel(fitTypes)
        histogram(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval)
    end
    ylabel('Count')
    xlabel('RMSE')
    % --- Quick ttest (sessions) --- 
    val1 = mdlFit.(fitTypes{1}).(cueLevel{iC}).(fitLevel).fval; % linear 
    val2 = mdlFit.(fitTypes{2}).(cueLevel{iC}).(fitLevel).fval; % linear2Hz
    varName = 'fval'; 
    [F.(varName).h, F.(varName).p, F.(varName).ci, F.(varName).stats] = ttest(val1,val2); 
    titleText = sprintf('Precue T%d: p = %0.4f, T = %0.2f', iC, F.(varName).p, F.(varName).stats.tstat); 
    title(titleText); 
end
legend(fitTypes)

% --- Save fig ---
figTitle = 'RMSEcompare'; 
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle)) 

%% Average sessions to subjects (if fit to sessions) 
% for iV = 1:2
%     for iC = 1:numel(cueLevel)
%         vals = [];
%         count = 1;
%         for iSession = 1:2:size(data.all.(fitLevel),2)
%             val = mean(fittedVal(iC,iSession:iSession+1));
%             vals(iC,count) = val;
%             count = count+1;
%         end
%         mdlFit.(cueLevel{iC}).(fitLevel).solutionSubject(iP,iV,iS)
%     end
% end

%% Polar plot of phase
iV = 4; iP = 1; 
% --- polar plot ---
figure
set(gcf,'Position',[100 100 400 400])
subplot 211
clear val 
for iC = 1:numel(cueLevel)
    for iS = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution,3)
        val(iC,iS) = squeeze( mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS) *2*pi/Fs ) ; 
        theta = [0 val(iC,iS)];
        rho = [0 1];
        polarplot(theta,rho,'Color',p.cueColors(iC,:))
        hold on
    end
end
subplot 212
for iS = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution,3)
    valDiff(iS) = val(1,iS) - val(2,iS);
    theta = [0 valDiff(iS)];
    rho = [0 1];
    polarplot(theta,rho,'Color',[0.5 0.5 0.5])
    hold on 
end
F(iV).p = circ_wwtest(val(1,:),val(2,:));

iV = 6; 
[F(iV).h, F(iV).p , F(iV).CI , F(iV).STATS] = ttest(valDiff); % 

% also try
polarhistogram(valDiff) % add 

%% Quick Ttests 
clear F 
for iV = 1:nVars
    for iF = 1:numel(fitTypes)
        clear val1 val2
        iP = 1;
        val1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(iP,iV,:);
        val2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(iP,iV,:);
        [F(iV).(fitTypes{iF}).h, F(iV).(fitTypes{iF}).p, F(iV).(fitTypes{iF}).ci, F(iV).(fitTypes{iF}).stats] = ttest(val1,val2);
    end
end

%% TTest on t(0) 
iV = 7; 
for iF = 1:numel(fitTypes)
    clear val1 val2
    val1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).yhat(:,1);
    val2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).yhat(:,1);
    [F(iV).(fitTypes{iF}).h, F(iV).(fitTypes{iF}).p, F(iV).(fitTypes{iF}).ci, F(iV).(fitTypes{iF}).stats] = ttest(val1,val2);
end

%% Prepare for export 
% Long format data for R
exportCSV = 1; 
for iF = 1:numel(fitTypes)
    clear V
    % fields = fieldnames(fittedP.(fitTypes{iF})); 
    switch fitTypes{iF}
        case 'linear2Hz'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope',...
                'amplitude','phase'};
        case 'linear'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope'};
    end
    count = 1;
    for iS = 1:10 % subjects
        for iC = 1:numel(cueLevel)
            for iSession = 1:2
                for iV = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution,2)
                    V.subject(count) = iS;
                    V.session(count) = iSession;
                    V.precue(count)  = iC;

                    % V.(tableHeaders{iV+3})(count) = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,iV);
                    V.(tableHeaders{iV+3})(count) = fittedP.(fitTypes{iF}).(tableHeaders{iV+3})(iS,iSession,iC); 
                    % V.(tableHeaders{iV+3})(count) = fittedP.(fitTypes{iF}).(fields{iV})(iS,iSession,iC); 
                end
                count = count+1;
            end
        end
    end

    % --- Export to csv for R ANOVA ---
    if exportCSV
        sz = [numel(V.(tableHeaders{1})) numel(tableHeaders)];

        varTypes = repmat("double",[1 numel(tableHeaders)]);
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

        for iT = 1:numel(tableHeaders)
            T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
        end

        % csvDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/model/ITPCfit1_constrainedPhase';
        csvName = sprintf('TANoise_ITPCFit_%s_%s',fitTypes{iF},dateStr);
        csvPath = sprintf('%s/%s.csv', figDir, csvName);
        writetable(T,csvPath)
    end
end

%% Prepare for export (OLD)
% Long format data for R
exportCSV = 1; 
for iF = 1:numel(fitTypes)
    clear V
    fields = fieldnames(fittedP.(fitTypes{iF})); 
    switch fitTypes{iF}
        case 'linear2Hz'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope',...
                'amplitude','phase','phaserad','freq',...
                'baseline_data','baseline_yhat'};
        case 'linear'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope',...
                'baseline_data','baseline_yhat'};
    end
    count = 1;
    for iS = 1:10 % subjects
        for iC = 1:numel(cueLevel)
            for iSession = 1:2
                for iV = 1:numel(fields)
                    V.subject(count) = iS;
                    V.session(count) = iSession;
                    V.precue(count)  = iC;

                    % V.(tableHeaders{iV+3})(count) = mdlFit.(cueLevel{iC}).(fitLevel).solution(iP,iV,iS);
                    V.(tableHeaders{iV+3})(count) = fittedP.(fitTypes{iF}).(fields{iV})(iS,iSession,iC); 
                end
                count = count+1;
            end
        end
    end

    % --- Export to csv for R ANOVA ---
    if exportCSV
        sz = [numel(V.(tableHeaders{1})) numel(tableHeaders)];

        varTypes = repmat("double",[1 numel(tableHeaders)]);
        T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders);

        for iT = 1:numel(tableHeaders)
            T.(tableHeaders{iT}) = V.(tableHeaders{iT})';
        end

        csvDir = '/Users/kantian/Dropbox/github/ta-meg-analysis2/model/ITPCfit1_constrainedPhase';
        csvName = sprintf('TANoise_ITPCFit_%s',fitTypes{iF});
        csvPath = sprintf('%s/%s.csv', csvDir, csvName);
        writetable(T,csvPath)
    end
end







