%% ITPC temporal expectation fit by precue (simultaneous fit of precue T1 and T2) 
% fmincon requires the optimization toolbox 
% set path to high level ta-meg-analysis2
% August 2023 Karen 

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian karen 

addpath(genpath(pwd))

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPC_yoked2/%s',user,dateStr); 
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
% === SETTINGS ===
fitLevel = 'session'; % 'session' 'group'
cueLevel = {'cueT1','cueT2'};
fitTypes = {'linear','linear2Hz'}; % 'linear2Hz'
o_bads = 1; % if 1, uses bads, otherwise uses fmincon
randCoeffs = 1; % uses random starting coefficients from the search space, otherwise uses hardcoded intial coefs
saveMdl = 1; % save mat of fitted parameters 
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
lb(idx) = -2; 
ub(idx) = 2; % 0.001; % arbitrarily large, for the scaling? 

% --- Amplitude --- 
idx = [find(contains(paramNames,'amplitude1')) find(contains(paramNames,'amplitude2'))]; 
lb(idx) = 0; 
ub(idx) = 1; 

% --- Phase shift ---
idx = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))]; 
lb(idx) = 0; 
ub(idx) = 2*pi; % 500 
if o_bads
    optionsBads = bads('defaults');
    optionsBads.PeriodicVars = idx;
    lb(idx) = -pi; 
    ub(idx) = pi;
end

% --- Frequency --- (fix to 2 Hz?) 
idx = find(contains(paramNames,'freq')); 
lb(idx) = 0.1; % *pi/(500); % 1.5 Hz 1.7
ub(idx) = 5; % *pi/(500); % 0.0016; % 2.5 Hz  2.3

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
            error('Specify starting coefficients manually')
        end

        clear dataFit1 dataFit2
        dataFit1 = data.(cueLevel{1}).(fitLevel)(tIdx,iS)';
        dataFit2 = data.(cueLevel{2}).(fitLevel)(tIdx,iS)';

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
    end
    if saveMdl
        % Save workspace and model fit
        filename = sprintf('ModelFit_YokedPrecueT1T2_%s',dateStr);
        save(filename,'-v7.3')
    end
end

%% --- Figure (data + fit) ---
for iF = 1:numel(fitTypes)
    for iS = 1:size(data.all.(fitLevel),2) % subjects
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
        clear dataFit1 dataFit2
        dataFit1 = data.(cueLevel{1}).(fitLevel)(:,iS)';
        dataFit2 = data.(cueLevel{2}).(fitLevel)(:,iS)';

        % --- Plot data ---
        plot(p.t,dataFit1,'LineWidth',1,'Color',p.cueColors(1,:))
        plot(p.t,dataFit2,'LineWidth',1,'Color',p.cueColors(2,:))

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
        [~,yhat1,yhat2] = meg_objectiveFunction2(fittedX,dataFit1(t),dataFit2(t),t,Fs,paramNames,fitTypes{iF});

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
            dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
            figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_%s_%s',sessionNames{iS},fitTypes{iF},dateStr);
            saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
        end
    end

    % lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast');
    % lgd.FontSize = 10;
end

%% Polar histogram (phase summary subject permutations
load('ModelFit_YokedPrecueT1T2_230930.mat')
mdlFitType = 'yoked'; 
meg_plotPhaseHistogram(mdlFit,paramNames,figDir,mdlFitType)

%% Polar histogram (phase summary subject permutations) 
edges = 0:pi/24:2*pi; 
for iS = 1:size(mdlFit.linear2Hz.(fitLevel).fval,2) % :size(data.all.(fitLevel),2) % subjects
    figure
    rhoPerm = 30; 
    % --- Perm w min fval --- 
    [minVal,idx] = min(  mdlFit.linear2Hz.(fitLevel).fval(:,iS)  );
    % --- Freq ---
    freq = mdlFit.linear2Hz.(fitLevel).solution(idx,iS,find(contains(paramNames,'freq'))); 
    % --- Phase in t ---
    valT_precueT1 = mdlFit.linear2Hz.(fitLevel).solution(:,iS,find(contains(paramNames,'phase1')));
    valT_precueT2 = mdlFit.linear2Hz.(fitLevel).solution(:,iS,find(contains(paramNames,'phase2')));
    % --- Calculate phase difference ---
    valR_precueT1 = t2rad( valT_precueT1, freq ); 
    valR_precueT2 = t2rad( valT_precueT2, freq ); 
    % --- Calculate phase difference ---
    valR_diff = circ_dist(valR_precueT1,valR_precueT2);

    % --- Plot fitted phases (radians) Precue T1 ---
    subplot 131
    rho = 1;
    polarhistogram(valR_precueT1,edges)
    hold on
    txt = sprintf('%0.2f deg, %0.2f\\pi rad', valT_precueT1(idx), valR_precueT1(idx)/pi);
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
    txt = sprintf('%0.2f deg, %0.2f\\pi rad', valT_precueT2(idx), valR_precueT2(idx)/pi);
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
    figTitle = sprintf('%s_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Phase_PermutationsHistogram_%d',...
        sessionNames{iS},...
        size(mdlFit.linear2Hz.(fitLevel).fval,1));
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Histogram of fitted freq 
for iS = 1:size(mdlFit.linear2Hz.(fitLevel).fval,2) % :size(data.all.(fitLevel),2) % subjects
    % --- Perm w min fval --- 
    [minVal,idx] = min(  mdlFit.linear2Hz.(fitLevel).fval(:,iS)  );
    freq(iS) = mdlFit.linear2Hz.(fitLevel).solution(idx,iS,find(contains(paramNames,'freq')));
end

figure
hold on
meg_figureStyle 
edges = 1:0.1:5; 
histogram(freq,edges)
ylabel('Count')
xlabel('Fitted frequency (Hz)')

% --- Save fig ---
figTitle = sprintf('Group_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Freq');
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% --- Figure (data + fit) (OLD) ---
% for iF = 1:numel(fitTypes)
%     figure
%     switch fitLevel
%         case 'session'
%             set(gcf,'Position',[100 100 1200 800])
%         case 'subject'
%             set(gcf,'Position',[100 100 800 800])
%         otherwise
%             error('Specify session- or subject-level fit')
%     end
% 
%     for iS = 1:size(data.all.(fitLevel),2) %
%         switch fitLevel
%             case 'session'
%                 subplot (size(data.all.(fitLevel),2)/4,4,iS)
%             case 'subject'
%                 subplot (size(data.all.(fitLevel),2)/2,2,iS)
%             otherwise
%                 error('Specify session- or subject-level fit')
%         end
%         hold on
% 
%         % --- Data (by precue) ---
%         dataFit1 = data.(cueLevel{1}).(fitLevel);
%         dataFit1 = dataFit1(:,iS)';
% 
%         dataFit2 = data.(cueLevel{2}).(fitLevel);
%         dataFit2 = dataFit2(:,iS)';
% 
%         % --- Plot data ---
%         plot(p.t, dataFit1,'LineWidth',2,'Color',p.cueColors(1,:))
%         plot(p.t, dataFit2,'LineWidth',2,'Color',p.cueColors(2,:))
% 
%         % --- Plot model fit ---
%         [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(fitLevel).fval(:,iS)  );
%         fittedX = mdlFit.(fitTypes{iF}).(fitLevel).solution(idx,:,iS);
%         switch fitTypes{iF}
%             case 'linear2Hz'
%                 paramNames = {'intercept1','slope1','amplitude1','phase1',...
%                     'intercept2','slope2','amplitude2','phase2',...
%                     'freq'};
%             case 'linear'
%                 paramNames = {'intercept1','slope1',...
%                     'intercept2','slope2'};
%         end
% 
%         [E,yhat1(iS,:),yhat2(iS,:)] = meg_objectiveFunction2(fittedX,dataFit1(toi),dataFit2(toi),t,Fs,paramNames,fitTypes{iF});
%         mdlFit.(fitTypes{iF}).(fitLevel).yhat1(iS,:) = yhat1(iS,:);
%         mdlFit.(fitTypes{iF}).(fitLevel).yhat2(iS,:) = yhat2(iS,:);
%         mdlFit.(fitTypes{iF}).(fitLevel).E(iS) = E; 
% 
%         fitColors = {'b','r'};
%         plot(p.t(toi), yhat1(iS,:),'LineWidth',1.5,'Color',fitColors{1})
%         plot(p.t(toi), yhat2(iS,:),'LineWidth',1.5,'Color',fitColors{2})
% 
%         % --- Plot event lines ---
%         for i = 1:numel(p.eventTimes)
%             xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
%         end
% 
%         % --- Format ---
%         meg_figureStyle
%         xlim([-100 2400])
%         if iS==size(data.all.(fitLevel),2)
%             xlabel('Time (ms)')
%             ylabel('ITPC')
%         end
%         % --- Titles --
%         switch fitLevel
%             case 'session'
%                 title(und2space(sessionNames{iS}))
%             case 'subject'
%                 title(und2space(subjectNames{iS}))
%             otherwise
%                 error('Specify session- or subject-level fit')
%         end
%         ax = gca;
%         ax.TitleFontSizeMultiplier = 0.8;
%         ax.TitleFontWeight = 'normal';
%     end
% 
%     % --- Save fig ---
%     dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
%     figTitle = sprintf('TANoise_ITPCFit_DataPrecueJoint_%s_%s',fitTypes{iF},dateStr);
%     saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
% end

%% --- Compare fvals (E) of yoked linear vs linear2Hz
clear val1 val2
figure
hold on
figureStyle
val1 = -1.*(mdlFit.linear2Hz.(fitLevel).fval); % full model 
val2 = -1.*(mdlFit.linear.(fitLevel).fval); % constrained model 
scatter(val1,val2,100,'filled','MarkerFaceColor',[0.5 0.5 0.5])
refline(1,0)
xlabel('-SSE Linear + 2Hz')
ylabel('-SSE Linear')
axis equal 
axis square 
dof = 4; % equal to the number of parameters that are constrained
clear stats 
[stats.h,stats.pValue,stats.stat,stats.cValue] = lratiotest(val1,val2,dof);

%% --- ttests on fitted model params (Precue T1 vs Precue T2) --- 
paramNamesAdd = paramNames; 
paramNamesAdd{end+1} = 'phaserad';
paramNamesAdd{end+1} = 't0'; 

clear fittedVal
for iF = 1:numel(fitTypes)
    for iV = 1:numel(paramNamesAdd)
        for iC = 1:numel(cueLevel)
            for iS = 1:size(data.all.(fitLevel),2)
                clear val
                if iV<=5
                    % Average to sessions ?
                    val = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS);
                    fittedVal(iV,iC,iS,iF) = val;
                end
            end
        end

        % Convert phase from t to rad, append
        switch paramNamesAdd{iV}
            case {'phase'} % circular stats
                fittedVal = deg2rad(fittedVal);
                val = (val*2*pi)/Fs; % covert from t to rad
                % fittedVal(iV,iC,iS,iF) = (fittedVal*2*pi)/Fs;
                % fittedVal(iV,iC,iS,iF) = val; 
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS) = fittedVal(iV,iC,iS,iF);
            % case {'intercept'}
            %     fittedX = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(1,:,iS);
            %     dataFit = data.(cueLevel{iC}).(fitLevel);
            %     dataFit = dataFit(:,iS);
            %     [E,yhat] = objectiveFunction(fittedX,dataFit,t,Fs);
            %     val = yhat(1); 
            %     fittedVal(iV,iC,iS,iF) = val; 
            %     mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS) = fittedVal(iV,iC,iS,iF);
        end
    end
end

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
for iV = 1:2
    for iC = 1:numel(cueLevel)
        vals = [];
        count = 1;
        for iSession = 1:2:size(data.all.(fitLevel),2)
            val = mean(fittedVal(iC,iSession:iSession+1));
            vals(iC,count) = val;
            count = count+1;
        end
        mdlFit.(cueLevel{iC}).(fitLevel).solutionSubject(iP,iV,iS)
    end
end

%% Figures of fitted params 
iV = 1;
vals = squeeze(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,:,:)); 
meg_plotFittedParam

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
clear V
tableHeaders = {'subject','session','precue','amplitude','intercept','slope','phaset','freq','phaserad'};
count = 1;
iPerm = 1;
for iS = 1:size(fittedVal(1,:),2)
    for iC = 1:numel(cueLevel)
        for iSession = 1:2
            for iV = 1:size(mdlFit.(cueLevel{iC}).(fitLevel).solution,2)
                V.subject(count) = iS;
                V.session(count) = iSession;
                V.precue(count)  = iC;

                V.(tableHeaders{iV+3})(count) = mdlFit.(cueLevel{iC}).(fitLevel).solution(iP,iV,iS);
                % disp(V.(tableHeaders{iV+3})(count))
            end
            count = count+1;
        end
    end
end

%% Export to csv for R ANOVA
sz = [numel(V.(tableHeaders{1})) numel(tableHeaders)]; 

varTypes = repmat("double",[1 numel(tableHeaders)]); 
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableHeaders); 

for iT = 1:numel(tableHeaders)
    T.(tableHeaders{iT}) = V.(tableHeaders{iT})'; 
end

csvDir = pwd; 
csvName = sprintf('TANoise_ITPCFit'); 
csvPath = sprintf('%s/%s.csv', csvDir, csvName); 
writetable(T,csvPath)







