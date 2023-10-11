% playground for testing 
% separate precue T1 & precue T2 fits 

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian (persona), karen (lab)

addpath(genpath(pwd))

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/model/ITPCsimulation_1/%s',user,dateStr); 
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
phase1 = 0; % change here 
% --- Data 2 params --- 
intercept2 = 0.3;
slope2 = 1; % 0.0001; 
amplitude2 = 0.1; 
phase2 = 1.5;
% --- Simulate data (in trial-relative time) --- 
clear dummyData
dummyData.cueT1 = ((slope1/1000 * p.t) + intercept1)  + ( amplitude1 * sin( (freq*pi/(Fs/2)) * (p.t + phase1 * 100 )) ); 
dummyData.cueT2 = ((slope2/1000 * p.t) + intercept2)  + ( amplitude2 * sin( (freq*pi/(Fs/2)) * (p.t + phase2 * 100 )) );

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
titleText = sprintf('Separate precue T1 and T2 fits\nintercept1 = %0.2f, slope1 = %0.5f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.5f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
    intercept1,slope1,...
    amplitude1,phase1,...
    intercept2,slope2,...
    amplitude2,phase2,...
    freq);
title(titleText)

ax = gca;
ax.TitleFontSizeMultiplier = 0.7;
ax.TitleFontWeight = 'normal';

% --- Save fig ---
figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_Separate_%s',dateStr);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% --- Fit settings --- 
paramNames = {'intercept','slope','amplitude','phase'};
nVars = numel(paramNames); 

% === SETTINGS === 
fitLevel = 'session'; % 'session' 'group' 
cueLevel = {'cueT1','cueT2'}; 
fitTypes = {'linear','linear2Hz'}; % 'linear2Hz'
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
lb(idx) = 0; 
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

% --- Nonlinear constraints --- 
nonlcon = []; 

% --- Options --- 
options = optimoptions('fmincon','Display','iter'); 

if o_bads
    optionsBads = bads('defaults');
    optionsBads.PeriodicVars = 4;
    mdlFit.lb(4) = -2.5;
    mdlFit.ub(4) = 2.5;
    mdlFit.hardlb = mdlFit.lb;
    mdlFit.hardub = mdlFit.ub;
end

%% --- Define initial search array --- 
nGrain = 100; 
clear x0s
x0s = NaN([nVars,nGrain]);
paramNames = {'intercept','slope','amplitude','phase'};
for i = 1:nVars
    % if strcmp(paramNames{i},'slope')
    %     x0s(i,:) = linspace(lb(i),0,nGrain); % slope max to 0? 
    % else
        x0s(i,:) = linspace(lb(i),ub(i),nGrain); % initial coefficients for search
    % end
end

%% --- Fit model ---
clear mdlFitDummy x0 x0Perm
nPerms = 1; % 100
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
    end
    
    for iC = 1:numel(cueLevel)
        clear dataFit
        dataFit = dummyData.(cueLevel{iC})(tIdx);
       
        for iF = 1:numel(fitTypes)
            % --- Save fit settings ---
            mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).x0(iP,:) = x0; % save starting coefficients
            if o_bads
                algo = 'bads';
            else
                algo = 'fmincon';
            end
            mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).algo = algo; % save optimization algorithm

            switch fitTypes{iF}
                case 'linear2Hz'
                    disp('Fitting linear + 2Hz ...')
                    paramNames = {'intercept','slope','amplitude','phase'};
                    % --- Define objective function ---
                    clear fun
                    fun = @(x)meg_objectiveFunction1(x,dataFit,t,Fs,paramNames,fitTypes{iF});
                    % --- Do fit ---
                    if o_bads
                        optionsBads = bads('defaults');
                        optionsBads.PeriodicVars = 4;
                        mdlFit.lb(4) = -2.5; 
                        mdlFit.ub(4) = 2.5; 
                        mdlFit.hardlb = mdlFit.lb; 
                        mdlFit.hardub = mdlFit.ub; 
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
            % --- Save fitted parameters --- 
            mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,:) = solution;
            % --- Save error --- 
            mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(iP) = fval;
            % --- Save predicted y --- 
            [E,yhat] = meg_objectiveFunction1(solution,dataFit,t,Fs,paramNames,fitTypes{iF});
            mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(iP,:) = yhat;
        end
    end
end

%% --- Figure (data + fit) ---
for iF = 1:numel(fitTypes)
    figure
    switch fitLevel
        case 'session'
            set(gcf,'Position',[100 100 450 300])
        case 'subject'
            set(gcf,'Position',[100 100 450 300])
        otherwise
            error('Specify session- or subject-level fit')
    end

    switch fitLevel
        case 'session'
            subplot (1,1,1)
        case 'subject'
            subplot (1,1,1)
        otherwise
            error('Specify session- or subject-level fit')
    end
    hold on
    
    for iC = 1:numel(cueLevel)
        % --- Data (by precue) ---
        dataFit = dummyData.(cueLevel{iC});

        % --- Plot data ---
        plot(p.t,dataFit,'LineWidth',1,'Color',p.cueColors(iC,:))

        % --- Plot model fit ---
        [minVal,idx] = min(  mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:)  );
        fittedX = mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,:);
        switch fitTypes{iF}
            case 'linear2Hz'
                paramNames = {'intercept','slope','amplitude','phase'};
            case 'linear'
                paramNames = {'intercept','slope'};
        end
        fitColors = {'b','r'};
        plot(t, mdlFitDummy.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(idx,:),'--','LineWidth',1,'Color',fitColors{iC})
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
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(2),...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(3),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(4),...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(2),...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(3),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(4),...
                2);
        case 'linear'
            titleText_Fitted = sprintf('Fitted: intercept1 = %0.2f, slope1 = %0.3f\nintercept2 = %0.2f, slope2 = %0.3f',...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(2),...
                mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(1),mdlFitDummy.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(2));
    end
    titleText1 = sprintf('Separate precue T1 and T2 fits (%s), fval = %0.2e',...
        mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
        mdlFitDummy.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval); 
    titleText_True = sprintf('True: intercept1 = %0.2f, slope1 = %0.3f, amplitude1 = %0.2f, phase1 = %0.2f\nintercept2 = %0.2f, slope2 = %0.3f, amplitude2 = %0.2f, phase2 = %0.2f\nfreq = %0.2f',...
        intercept1, slope1, amplitude1, phase1,...
        intercept2, slope2, amplitude2, phase2,...
        freq); 
    titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted,titleText_True); 
    title(titleText)
    
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.6;
    ax.TitleFontWeight = 'normal';

    lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast'); 
    lgd.FontSize = 10; 

    % --- Save fig ---
    dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
    figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_Separate_%s_%s',fitTypes{iF},dateStr);
    saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))
end

%% Polar plot of fitted phase
figure
subplot 121 
meg_figureStyle
% --- Calculate phase difference --- 
% valT_diff = abs(mdlFitDummy.linear2Hz.(cueLevel{1}).(fitLevel).solution(4) - mdlFitDummy.linear2Hz.(cueLevel{2}).(fitLevel).solution(4));
valR_diff = abs(circ_dist(mdlFitDummy.linear2Hz.(cueLevel{1}).(fitLevel).solution(4), mdlFitDummy.linear2Hz.(cueLevel{2}).(fitLevel).solution(4))); 
valT_diff = rad2t(valT_diff);
% --- Plot fitted phases (radians) --- 
for iC = 1:numel(cueLevel)
    rho = 1; 
    valRad = mdlFitDummy.linear2Hz.(cueLevel{iC}).(fitLevel).solution(4);
    valT = rad2t(valT); 
    polarscatter(valRad,rho,100,'filled','MarkerFaceColor',p.cueColors(iC,:))
    hold on 
    polarplot([0 valRad],[0 rho],'Color',p.cueColors(iC,:))
    txt = sprintf('%0.2f deg = %0.2f rad',valT,valRad); 
end
% --- Plot phase difference --- 
polarplot([0 valR_diff],[0 rho],'Color',[0.5 0.5 0.5])
% --- Formatting --- 
pax = gca; 
pax.ThetaAxisUnits = 'radians';
title('Fitted phase (rad)')

subplot 122 
meg_figureStyle
% --- Plot fitted phases (t space) --- 
for iC = 1:numel(cueLevel)
    rho = 1; 
    valRad = mdlFitDummy.linear2Hz.(cueLevel{iC}).(fitLevel).solution(4);
    valT = rad2t(valT); 
    polarscatter(valRad,rho,100,'filled','MarkerFaceColor',p.cueColors(iC,:))
    hold on 
    polarplot([0 valRad],[0 rho],'Color',p.cueColors(iC,:))
    txt = sprintf('%0.2f deg = %0.2f rad',valT,valRad); 
end
% --- Plot phase difference --- 
polarplot([0 valR_diff],[0 rho],'Color',[0.5 0.5 0.5])
% --- Formatting --- 
pax = gca; 
pax.ThetaAxisUnits = 'radians';
ticks = 0:pi/6:2*pi; 
thetaticks(ticks);
thetaticklabels(rad2t(ticks))
title('Fitted phase (t/100)')

% --- Save fig ---
dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
figTitle = sprintf('Dummy_TANoise_ITPCFit_DataPrecue_Separate_Linear2Hz_Phase_%s',dateStr);
saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% Check conversion of phase into polar plot
val = -2.5:0.5:2.5;  % fittedP.(fitTypes{iF}).(paramNames{iV});
% val = (val*2*pi)/Fs; % if -500 to 500 bounds
% val = (val*100)*(4*freq*pi)/Fs;
% test conversion
% val = t2rad(val); 
% 
% figure
% rho = ones(size(val));
% polarscatter(val,rho,50,'filled','MarkerFaceColor','k')
% hold on
% 
% pax = gca; 
% pax.ThetaAxisUnits = 'radians';
% 
% dateStr = datetime('now','TimeZone','local','Format','yyMMdd_hhmm');
% figTitle = sprintf('Dummy_TANoise_phaseTest_%s',dateStr);
% saveas(gcf,sprintf('%s/%s.png', figDir, figTitle))

%% Histogram of data? 
figure 
meg_figureStyle
hold on
histogram(dummyData.cueT1)
% for i = 1:size(dummyData.cueT1,2)
%     plot(dummyData.cueT1(i))
% end




