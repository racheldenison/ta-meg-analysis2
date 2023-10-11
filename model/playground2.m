% playground for testing 
% yoked precue T1 & precue T2 fits 

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian karen 

addpath(genpath(pwd))

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPC2_playground/%s',user,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% --- Load data parameters ---
p = meg_params('TANoise_ITPCsession8');

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

% --- Data settings ---
foi = 20; % frequency of interest, Hz
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 
Fs = 1000; % sampling frequency 
fitLevel = 'session'; 

% -- MEG settings --- 
expt = 'TANoise'; 
[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

% --- Extract data into precue conds at desired foi --- 
cueNames = {'all','cueT1','cueT2'};
fitLevels = {'session','subject','group'}; 

%% --- Generate simulated data ---
Fs = 1000; 
freq = 2; 
% --- Data 1 params --- 
intercept1 = 0.3; 
slope1 = 1; % 0.0001; % 0.0001 more reflects actual data 
amplitude1 = 0.1; 
phase1_t = rad2t(1*pi,freq,Fs); 
phase1_rad = t2rad(phase1_t,freq,Fs); % change here, try editing to rad space 
% --- Data 2 params --- 
intercept2 = 0.3;
slope2 = 1.01; % 0.0001; 
amplitude2 = 0.2; 
phase2_t = rad2t(1.8*pi,freq,Fs); % phase is specified in t/100 for parameter scaling 
phase2_rad = t2rad(phase2_t,freq,Fs); 
% --- Simulate data (in trial-relative time) --- 
clear dummyData
% dummyData.cueT1 = ((slope1/1000 * p.t) + intercept1)  + ( amplitude1 * sin( (freq*pi/(Fs/2)) * (p.t + phase1 * 100 )) ); 
% dummyData.cueT2 = ((slope2/1000 * p.t) + intercept2)  + ( amplitude2 * sin( (freq*pi/(Fs/2)) * (p.t + phase2 * 100 )) );
dummyData.cueT1 = ((slope1/1000 * p.t) + intercept1)  + ( amplitude1 * sin( (freq*pi/(Fs/2)) * (p.t + phase1_t * 100 )) ); 
dummyData.cueT2 = ((slope2/1000 * p.t) + intercept2)  + ( amplitude2 * sin( (freq*pi/(Fs/2)) * (p.t + phase2_t * 100 )) );

% --- Add noise ---
addNoise = 0; 
if addNoise
    noise1 = (rand(size(dummyData.cueT1))-0.5) * 0.5;
    noise2 = (rand(size(dummyData.cueT2))-0.5) * 0.2;

    dummyData.cueT1 = dummyData.cueT1 + noise1;
    dummyData.cueT2 = dummyData.cueT2 + noise2;
end

%% Plot data 
figure
switch fitLevel
    case 'session'
        set(gcf,'Position',[100 100 500 300])
    case 'subject'
        set(gcf,'Position',[100 100 500 300])
    otherwise
        error('Specify session- or subject-level fit')
end
hold on 

% --- Plot data ---
plot(p.t,dummyData.cueT1,'LineWidth',2)
plot(p.t,dummyData.cueT2,'LineWidth',2)
xlabel('Time (ms)')
ylabel('Simulated ITPC')

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
titleText = sprintf('Yoked precue T1 and T2 fits\nintercept1 = %0.2f, slope1 = %0.5f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.5f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
    intercept1,slope1,...
    amplitude1,phase1_rad,...
    intercept2,slope2,...
    amplitude2,phase2_rad,...
    freq);
title(titleText)

ax = gca;
ax.TitleFontSizeMultiplier = 0.7;
ax.TitleFontWeight = 'normal';

% --- Save fig ---
figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_%s');
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% --- Fit settings --- 
% === SETTINGS ===
fitLevel = 'session'; % 'session' 'group'
cueLevel = {'cueT1','cueT2'};
fitTypes = {'linear','linear2Hz'}; % 'linear2Hz'
o_bads = 1; % if 1, uses bads, otherwise uses fmincon
randCoeffs = 1; % uses random starting coefficients from the search space, otherwise uses hardcoded intial coefs
saveMdl = 0; % save mat of fitted parameters 
% ================

% 1 for precue T1, 2 for precue T2
paramNames = {'intercept1','slope1','amplitude1','phase1',...
              'intercept2','slope2','amplitude2','phase2',...
              'freq'};
nVars = numel(paramNames); 

clear mdlFit lb ub 
mdlFit.A = []; % equality contraints
mdlFit.b = []; 
mdlFit.Aeq = []; % inequality contraints
mdlFit.beq = []; 

% === Define upper and lower bounds ===

% --- Intercept ---
idx = [find(contains(paramNames,'intercept1')) find(contains(paramNames,'intercept2'))]; 
lb(idx) = 0; 
ub(idx) = 1; 

% --- Slope --- 
idx = [find(contains(paramNames,'slope1')) find(contains(paramNames,'slope2'))]; 
lb(idx) = 0; 
ub(idx) = 2; % 0.001; % arbitrarily large, for the scaling? 

% --- Amplitude --- 
idx = [find(contains(paramNames,'amplitude1')) find(contains(paramNames,'amplitude2'))]; 
lb(idx) = 0; 
ub(idx) = 1; 

% --- Phase shift ---
idx = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))]; 
lb(idx) = 0; % -500 idx = find(contains(paramNames,'slope')); 
ub(idx) = 2*pi; % 5, 500 
if o_bads
    optionsBads = bads('defaults');
    optionsBads.PeriodicVars = idx;
    lb(idx) = -pi; % -2.5 
    ub(idx) = pi; % 2.5 
end

% --- Frequency --- (fix to 2 Hz?) 
idx = find(contains(paramNames,'freq')); 
lb(idx) = 1.9; % *pi/(500); % 1.5 Hz 1.7
ub(idx) = 2.1; % *pi/(500); % 0.0016; % 2.5 Hz  2.3

mdlFit.lb = lb; 
mdlFit.ub = ub; 

% --- Nonlinear constraints --- 
nonlcon = []; 

% --- Options --- 
options = optimoptions('fmincon','Display','iter'); 

if o_bads
    mdlFit.hardlb = mdlFit.lb;
    mdlFit.hardub = mdlFit.ub;
end

%% --- Define initial search array --- 
nGrain = 100; 
clear x0s
x0s = NaN([nVars,nGrain]);
for i = 1:nVars
    x0s(i,:) = linspace(mdlFit.lb(i),mdlFit.ub(i),nGrain); % initial coefficients for search
end

%% --- Fit model ---
clear x0 x0Perm
nPerms = 10; % 10 100
for iS = 1 % :size(data.all.(fitLevel),2) % subjects
    for iP = 1:nPerms
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
            error('Specify starting coefficients manually')
        end

        clear dataFit1 dataFit2
        dataFit1 = dummyData.(cueLevel{1})(tIdx);
        dataFit2 = dummyData.(cueLevel{2})(tIdx);

        for iF = 1:numel(fitTypes)
            % --- Save fit settings ---
            mdlFit.(fitTypes{iF}).(fitLevel).x0(iP,iS,:) = x0; % save starting coefficients
            if o_bads
                algo = 'bads';
            else
                algo = 'fmincon';
            end
            mdlFit.(fitTypes{iF}).(fitLevel).algo = algo; % save optimization algorithm
            
            paramNames = {'intercept1','slope1','amplitude1','phase1',...
                'intercept2','slope2','amplitude2','phase2','freq'};
            switch fitTypes{iF}
                case 'linear2Hz'
                    disp('Fitting linear + 2Hz ...')
                    % --- Define objective function ---
                    clear fun
                    fun = @(x)meg_objectiveFunction2(x,dataFit1,dataFit2,t,Fs,paramNames,fitTypes{iF});
                    % --- Do fit ---
                    if o_bads
                        [solution,fval,exitflag,output] = bads(fun, x0, mdlFit.hardlb, mdlFit.hardub, mdlFit.lb, mdlFit.ub,[],optionsBads);
                    else
                        [solution,fval,exitflag,output] = fmincon(fun, x0, mdlFit.A, mdlFit.b, mdlFit.Aeq, mdlFit.beq, mdlFit.lb, mdlFit.ub);
                    end
                case 'linear'
                    disp('Fitting linear ...')
                    idx = [find(contains(paramNames,'intercept1')) find(contains(paramNames,'slope1')) ...
                        find(contains(paramNames,'intercept2')) find(contains(paramNames,'slope2'))]; 
                    paramNamesLinear = {'intercept1','slope1',...
                        'intercept2','slope2'};
                    % --- Define objective function ---
                    clear fun
                    fun = @(x)meg_objectiveFunction2(x,dataFit1,dataFit2,t,Fs,paramNamesLinear,fitTypes{iF});
                    % --- Do fit ---
                    if o_bads
                        [solution,fval,exitflag,output] = bads(fun, x0(idx), mdlFit.hardlb(idx), mdlFit.hardub(idx), mdlFit.lb(idx), mdlFit.ub(idx));
                    else
                        [solution,fval,exitflag,output] = fmincon(fun, x0(idx), mdlFit.A, mdlFit.b, mdlFit.Aeq, mdlFit.beq, mdlFit.lb(idx), mdlFit.ub(idx));
                    end
            end
            % --- Save fitted parameters ---
            mdlFit.(fitTypes{iF}).(fitLevel).solution(iP,iS,:) = solution;
            % --- Save error ---
            mdlFit.(fitTypes{iF}).(fitLevel).fval(iP,iS) = fval;
            % --- Save predicted y ---
            % clear E yhat
            % [E,yhat] = meg_objectiveFunction1(solution,dataFit,t,Fs,paramNames,fitTypes{iF});
            % mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(iP,iS,:) = yhat;
        end
        if saveMdl
            % Save workspace and model fit
            filename = sprintf('ModelFit_YokedPrecueT1T2_%s',dateStr);
            save(filename,'-v7.3')
        end
    end
end

%% --- Figure (data + fit) ---
for iF = 1:numel(fitTypes)
    for iS = 1 % :size(data.all.(fitLevel),2) % subjects
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
        % --- Plot data ---
        plot(p.t,dummyData.(cueLevel{iC}),'LineWidth',1,'Color',p.cueColors(iC,:))
    end

    % --- Find model fit with lowest fval ---
    [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(fitLevel).fval(:,iS)  );
    fittedX = squeeze(mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,:));
    % --- Plot model fit ---
                switch fitTypes{iF}
                    case 'linear2Hz'
                        paramNames = {'intercept1','slope1','amplitude1','phase1',...
                            'intercept2','slope2','amplitude2','phase2','freq'};
                    case 'linear'
                        paramNames= {'intercept1','slope1',...
                            'intercept2','slope2'};
                end
    [~,yhat1,yhat2] = meg_objectiveFunction2(fittedX,dataFit1,dataFit2,t,Fs,paramNames,fitTypes{iF}); 

    fitColors = {'b','r'};
    plot(t,yhat1,'--','LineWidth',1,'Color',fitColors{1})
    plot(t,yhat2,'--','LineWidth',1,'Color',fitColors{2})

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
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,1),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,2),...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,3),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,4),...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,5),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,6),...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,7),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,8),...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,9));
        case 'linear'
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f\nintercept2 = %0.2f, slope2 = %0.3f',...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,1),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,2),...
                mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,3),mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,iS,4));
    end
    titleText1 = sprintf('Yoked precue T1 and T2 fits (%s), fval = %0.2e, %s',...
        mdlFit.(fitTypes{iF}).(fitLevel).algo,...
        mdlFit.(fitTypes{iF}).(fitLevel).fval(iS),...
        und2space(sessionNames{iS})); 
    titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted); 
    title(titleText)
    
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.6;
    ax.TitleFontWeight = 'normal';

    % --- Save fig ---
    if saveFigs
        dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
        figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_%s_%s',sessionNames{iS},fitTypes{iF},dateStr);
        saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    end
    end

    % lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast');
    % lgd.FontSize = 10;
end

%% Polar histogram (phase summary subject permutations) 
saveFigs = 0; 
edges = 0:pi/24:2*pi; 
for iS = 1:size(mdlFit.linear2Hz.(fitLevel).fval,2) % :size(data.all.(fitLevel),2) % subjects
    figure
    rhoPerm = 100; 
    % --- Perm w min fval --- 
    [minVal,idx] = min(  mdlFit.linear2Hz.(fitLevel).fval(:,iS)  );
    % --- Freq ---
    freq = mdlFit.linear2Hz.(fitLevel).solution(idx,iS,find(contains(paramNames,'freq'))); 
    % --- Phase in t ---
    valR_precueT1 = mdlFit.linear2Hz.(fitLevel).solution(:,iS,find(contains(paramNames,'phase1')));
    valR_precueT2 = mdlFit.linear2Hz.(fitLevel).solution(:,iS,find(contains(paramNames,'phase2')));
    % --- Calculate phase difference ---
    valT_precueT1 = rad2t( valR_precueT1, freq ); 
    valT_precueT2 = rad2t( valR_precueT2, freq ); 
    % --- Calculate phase difference ---
    valR_diff = circ_dist(valR_precueT1,valR_precueT2);

    % --- Plot fitted phases (radians) Precue T1 ---
    subplot 131
    rho = 1;
    polarhistogram(valR_precueT1,edges)
    hold on
    txt = sprintf('%0.2f ms/100, %0.2f\\pi rad', valT_precueT1(idx), valR_precueT1(idx)/pi);
    % --- Plot fitted phase w min fval - Precue T1 ---
    polarplot([0 valR_precueT1(idx)],[0 rhoPerm],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(sprintf('Precue T1\n%s',txt))

    % --- Plot fitted phases (radians) Precue T2 ---
    subplot 132
    rho = 1;
    polarhistogram(valR_precueT2,edges)
    hold on
    txt = sprintf('%0.2f ms/100, %0.2f\\pi rad', valT_precueT2(idx), valR_precueT2(idx)/pi);
    % --- Plot fitted phase w min fval - Precue T2 ---
    polarplot([0 valR_precueT2(idx)],[0 rhoPerm],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    title(sprintf('Precue T2\n%s',txt))

    % --- Plot phase difference ---
    subplot 133
    polarhistogram(valR_diff,edges)
    hold on
    % --- Plot fitted phase w min fval - Precue T1-Precue T2 ---
    polarplot([0 valR_diff(idx)],[0 rhoPerm],'r')
    % --- Formatting ---
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    txt = sprintf('%0.2f\\pi rad', valR_diff(idx)/pi);
    title(sprintf('Precue T1 - Precue T2\n%s',txt))

    sgtitle(sprintf('Fitted phase %s, nPerm = %d',...
        und2space(sessionNames{iS}),...
        size(mdlFit.linear2Hz.(fitLevel).fval,1)))

    % --- Save fig ---
    if saveFigs
        figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Phase_PermutationsHistogram_%d',...
            sessionNames{iS},...
            size(mdlFit.linear2Hz.(fitLevel).fval,1));
        saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
    end
end


%% --- Figure (data + fit) ---
% for iF = 1:numel(fitTypes)
%     figure
%     switch fitLevel
%         case 'session'
%             set(gcf,'Position',[100 100 500 300])
%         case 'subject'
%             set(gcf,'Position',[100 100 500 300])
%         otherwise
%             error('Specify session- or subject-level fit')
%     end
% 
%     switch fitLevel
%         case 'session'
%             subplot (1,1,1)
%         case 'subject'
%             subplot (1,1,1)
%         otherwise
%             error('Specify session- or subject-level fit')
%     end
%     hold on
% 
%     % --- Data (by precue) ---
%     dataFit1 = dummyData1;
%     dataFit2 = dummyData2;
% 
%     % --- Plot data ---
%     plot(p.t(toi), dataFit1,'LineWidth',1,'Color',p.cueColors(1,:))
%     plot(p.t(toi), dataFit2,'LineWidth',1,'Color',p.cueColors(2,:))
% 
%     % --- Plot model fit ---
%     [minVal,idx] = min(  mdlFitDummy.(fitTypes{iF}).(fitLevel).fval(:)  );
%     fittedX = mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(idx,:);
%     switch fitTypes{iF}
%         case 'linear2Hz'
%             paramNames = {'intercept1','slope1','amplitude1','phase1',...
%                 'intercept2','slope2','amplitude2','phase2',...
%                 'freq'};
%         case 'linear'
%             paramNames = {'intercept1','slope1',...
%                 'intercept2','slope2'};
%     end
% 
%     [E,yhat1,yhat2] = meg_objectiveFunction2(fittedX,dataFit1,dataFit2,t,Fs,paramNames,fitTypes{iF});
%     mdlFitDummy.(fitTypes{iF}).(fitLevel).yhat1 = yhat1;
%     mdlFitDummy.(fitTypes{iF}).(fitLevel).yhat2 = yhat2;
%     mdlFitDummy.(fitTypes{iF}).(fitLevel).E = E;
% 
%     fitColors = {'b','r'};
%     plot(p.t(toi), yhat1,':','LineWidth',3,'Color',fitColors{1})
%     plot(p.t(toi), yhat2,':','LineWidth',3,'Color',fitColors{2})
% 
%     % --- Plot event lines ---
%     for i = 1:numel(p.eventTimes)
%         xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
%     end
% 
%     % --- Format ---
%     meg_figureStyle
%     xlim([-100 2400])
%     xlabel('Time (ms)')
%     ylabel('ITPC')
% 
%     % --- Titles --
%     switch fitTypes{iF}
%         case 'linear2Hz'
%             titleText = sprintf('Yoked precue T1 and T2 fits\nintercept1 = %0.2f, slope1 = %0.5f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.5f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
%                 mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(2),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(3),...
%                 mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(4),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(5),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(6),...
%                 mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(7),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(8),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(9));
%         case 'linear'
%             titleText = sprintf('Yoked precue T1 and T2 fits\nintercept1 = %0.2f, slope1 = %0.5f\nintercept2 = %0.2f, slope2 = %0.5f',...
%                 mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(2),mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(3),...
%                 mdlFitDummy.(fitTypes{iF}).(fitLevel).solution(4));
%     end
%     title(titleText)
% 
%     ax = gca;
%     ax.TitleFontSizeMultiplier = 0.7;
%     ax.TitleFontWeight = 'normal';
% 
%     % --- Save fig ---
%     dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
%     figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_yoked_%s_%s',fitTypes{iF},dateStr);
%     saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
% end

%% --- Fit model (dummy test old) --- 
clear dataFit fun solution x0
nPerms = 1; % 100
for iP = 1:nPerms
    % --- Randomly pick starting coefficients from search grid ---
    for iV = 1:nVars
        vals = randperm(nGrain);
        idx(iV) = vals(1);
        x0(iV) = x0s(iV,idx(iV));
    end
    fun = @(x)playground_objectiveFunction(x,dummyData1,dummyData2,t,paramNames);
    % --- Fit ---
    [solution,fval,exitflag,output] = fmincon(fun, x0, mdlFit.A, mdlFit.b, mdlFit.Aeq, mdlFit.beq, mdlFit.lb, mdlFit.ub);
end

%% Plot data and model 
figure 
hold on 
figureStyle 

% --- Plot data ---
plot(t,dummyData1,'LineWidth',1,'Color',p.cueColors(1,:))
plot(t,dummyData2,'LineWidth',1,'Color',p.cueColors(2,:))

% --- Plot model fit ---
[E,yhat1,yhat2] = playground_objectiveFunction(solution,dummyData1,dummyData2,t,paramNames);
plot(t,yhat1,':','LineWidth',2,'Color','b')
plot(t,yhat2,':','LineWidth',2,'Color','r')












