% meg_fitAnticipatory1_freeFreq
% all trials (or single precue condition), free frequency 
% November 2023 Karen

%% Settings
saveFigs = 1; 
user = 'kantian'; % kantian (persona), karen (lab)

workingDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis-model/model_anticipatory',user); 
addpath(genpath(workingDir)) % ta-meg-analysis-model/model_anticipatory
addpath(genpath('../../ta-meg-analysis2')) 

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('%s/ITPCfit_1_freeFreq/%s',workingDir,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

%% Prepare data 
% --- Load ITPC data ---
load(sprintf('/Users/%s/Dropbox/Data/TANoise/MEG/Group/mat/groupA_ITPCspectrogram_byAtt.mat',user)) % A variable 
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
paramNames = {'intercept','slope','amplitude','phase','freq'};
nVars = numel(paramNames); 

% === SETTINGS ===
fitLevel = 'session'; % 'session' 'subject' 'group'
cueLevel = {'all','cueT1','cueT2'}; % 'cueT1','cueT2'
fitTypes = {'linear','linear2Hz'}; % 'linear2Hz' 'linear'
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
idx = find(contains(paramNames,'freq')); 
lb(idx) = 0.1; % *pi/(500); % 1.5 Hz 1.7
ub(idx) = 5; % *pi/(500); % 0.0016; % 2.5 Hz  2.3

mdlFit.lb = lb; 
mdlFit.ub = ub; 

mdlFit.hardlb = mdlFit.lb;
mdlFit.hardub = mdlFit.ub;

% --- Nonlinear constraints --- 
nonlcon = []; 

% --- Options --- 
options = optimoptions('fmincon','Display','iter'); 

%% --- Define initial search array --- 
nGrain = 100; 
clear x0s
x0s = NaN([nVars,nGrain]);
paramNames = {'intercept','slope','amplitude','phase','freq'};
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
                        paramNames = {'intercept','slope','amplitude','phase','freq'};
                        % --- Define objective function ---
                        clear fun
                        fun = @(x)meg_objectiveFunction1_freeFreq(x,dataFit,t,Fs,paramNames,fitTypes{iF});
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
    filename = sprintf('ModelFit_1_freeFreq_%s',dateStr);
    save(filename,'-v7.3')
end

%% Load data (if regenerating) 
load('/Users/kantian/Dropbox/github/ta-meg-analysis-model/model_anticipatory/ModelFit_1_freeFreq_GroupTS_231102.mat')

%% --- Figure (data + fit) FOR MANUSCRIPT !!! ---
titleVis = 0; % if title vis off, then will plot for appropriate manuscript size 
showN = 1; % show n = X annotation 
figFormat = 'svg'; % svg 
plotErrorBars = 1; % turn off for subject-level plots 
restrictYLim = 0; % turn on for group-level manuscript matching ylims 

for iF = 1:numel(fitTypes)
    for iS = 1:size(data.all.(fitLevel),2) % sessions
    figure
    switch fitLevel
        case {'session','subject','group'}
            set(gcf,'Position',[100 100 500 300]) %  [100 100 600 400]
        otherwise
            error('Specify session- or subject-level fit')
    end
    hold on

    % --- Format ---
    meg_figureStyle
    [style, colors] = meg_manuscriptStyle;
    xlim([-300 2400]); % xlim([-100 2400])
    if restrictYLim
        ylim([0.25 0.4])
        yticks(0.25:0.05:0.4)
    end
    xlabel('Time (ms)')
    ylabel('ITPC')

    % --- Plot event lines ---
    for i = 1:numel(p.eventTimes)
        xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
    end
    
    for iC = 1:numel(cueLevel)
        % --- Plot data error bars ---
        % Requires https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
        sampling = 1:10:7001;
        errorBarType = 'EB-TS-baselineNorm'; % 'EB-TS' 'EB-fit' 'EB-TS-baselineNorm'
        if plotErrorBars
            includeIdx = 1:10;
            groupVals = data.(cueLevel{iC}).group; 
            vals = data.(cueLevel{iC}).subject; 
            switch errorBarType
                % case 'EB-TS'
                %     t = p.t(sampling);
                %     line_allTrials = shadedErrorBar(t, groupVals(sampling), std(vals(sampling,includeIdx),[],2)/sqrt(numel(includeIdx)),...
                %         'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',1);
                case 'EB-TS-baselineNorm'
                    baselineTOI = -100:-1; % why is this -200? 
                    baselineTOIIdx = find(p.t==baselineTOI(1)):find(p.t==baselineTOI(end));
                    for i = 1:size(vals,2) % average per subject
                        baselineS(i) = mean(vals(baselineTOIIdx,i),1);
                        normVals(:,i) = vals(:,i)./baselineS(i);
                    end
                    % sem in units of percent change
                    sem = std(normVals(sampling,includeIdx),[],2)/sqrt(numel(includeIdx));
                    % scale sem back to ITPC units
                    baselineG = mean(groupVals(baselineTOIIdx));
                    semScaled = sem*baselineG;
                    % Plot
                    t_plot = p.t(sampling);
                    line_allTrials = shadedErrorBar(t_plot, groupVals(sampling), semScaled,...
                        'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',0);
                    line_allTrials.patch.FaceColor = colors.lightPurple; 
                    line_allTrials.patch.FaceAlpha = 1; 
                    for iL = 1:2
                        line_allTrials.edge(iL).Color = colors.mediumPurple; 
                        line_allTrials.edge(iL).LineWidth = 0.2; 
                    end
                % case 'EB-fit'
                %     t = p.t(toi);
                %     % Collect subject-level fits into single matrix
                %     for i = 1:numel(c_subs)
                %         sub_fits(:,i) = c_subs(i).y1_est;
                %     end
                %     line_allTrials = shadedErrorBar(t, y1_est, std(sub_fits(:,includeIdx),[],2)/sqrt(numel(includeIdx)),...
                %         'lineProps', {'MarkerFaceColor','k', 'LineWidth', 0.2, 'Color','k'}, 'transparent',1);
            end
        end

        % --- Data (by precue) ---
        dataFit = data.(cueLevel{iC}).(fitLevel)(:,iS)';

        % --- Plot data ---
        plot(p.t(sampling),dataFit(sampling),'LineWidth',2,'Color',colors.mediumPurple)

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
                paramNames = {'intercept','slope','amplitude','phase','freq'};
            case 'linear'
                paramNames = {'intercept','slope'};
        end
        [~,yhat] = meg_objectiveFunction1_freeFreq(fittedX,dataFit(toi),t,Fs,paramNames,fitTypes{iF});
        plot(t,yhat,'--','LineWidth',2,'Color',colors.darkPurple)
    end

    % --- Plot shaded baseline window ---
    xl = xlim;
    yl = ylim;
    x = [p.t(tIdx(1)) p.t(tIdx(end)) p.t(tIdx(end)) p.t(tIdx(1))];
    % y = [yl(1) yl(1) yl(2) yl(2)]; % [yl(1) yl(1) yl(2) yl(2)]
    % patch(x,y,lightgrey,'EdgeColor',lightgrey)
    yScale = yl(2)-yl(1); 
    ySize = (0.02*0.15)*(yScale/0.15);
    y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize]; % [yl(1) yl(1) yl(2) yl(2)]
    
    patch(x,y,colors.mediumgrey,'EdgeColor',colors.mediumgrey)
    % plot([x(1) x(2)],[y(1) y(1)],'Color',mediumgrey,'LineWidth',30)

    % --- Plot pre-precue baseline check --- 
    x = [p.t(tIdx(1))-100 p.t(tIdx(1)) p.t(tIdx(1)) p.t(tIdx(1))-100];
    y = [yl(1) yl(1) yl(1)+ySize yl(1)+ySize];
    patch(x,y,colors.lightgrey,'EdgeColor',colors.lightgrey)

    if showN
        % --- Add n annotation ---
        nStr = sprintf('n = %d x 2 sessions',size(data.all.subject,2));
        nStrTxt = text(0.945*xl(2),yl(1)+0.003,nStr,'HorizontalAlignment','right','VerticalAlignment','bottom');
        nStrTxt.FontSize = 14;
        nStrTxt.FontName = 'Helvetica-Light';
    end

    % set(gca,'children',flipud(get(gca,'children'))) % change order of plotted objects 
    set(gca ,'Layer', 'Top') % bring axes to front

    % --- Titles --
    if titleVis
        switch fitTypes{iF}
            case 'linear2Hz'
                titleText_Fitted = sprintf('Fitted: intercept = %0.2f, slope = %0.3f, amplitude = %0.2f, phase = %0.2f\\pi\nfreq = %0.2f',...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2),...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,3),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,4)/pi,...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,5) );
            case 'linear'
                titleText_Fitted = sprintf('Fitted: intercept = %0.2f, slope = %0.3f',...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,1),mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS,2));
        end
        switch fitLevel
            case {'sessions','subjects'}
                titleText1 = sprintf('Fits to all trials (%s), fval = %0.2e, %s',...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(iS),...
                    und2space(sessionNames{iS}));
            case 'group'
                titleText1 = sprintf('Fits to all trials (%s), fval = %0.2e, n = %d',...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).algo,...
                    mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).fval(iS),...
                    size(data.all.subject,2));
        end
        titleText = sprintf('%s\n%s\n%s',titleText1,titleText_Fitted);
        title(titleText)

        ax = gca;
        ax.TitleFontSizeMultiplier = 0.6;
        ax.TitleFontWeight = 'normal';
    end

    % --- Save fig ---
    dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
    switch fitLevel
        case 'session'
            figTitle = sprintf('%s_TANoise_ITPCFit_AllTrials_%s_%s',sessionNames{iS},fitTypes{iF},dateStr);
        case 'subject'
            figTitle = sprintf('%s_TANoise_ITPCFit_AllTrials_%s_%s',subjectNames{iS},fitTypes{iF},dateStr);
        case 'group'
            figTitle = sprintf('Group_TANoise_ITPCFit_AllTrials_%s_%s',fitTypes{iF},dateStr);
    end
    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
    end

    % lgd = legend({'Precue T1 data','Precue T1 fit','Precue T2 data','Precue T2 fit'},'Location','southeast');
    % lgd.FontSize = 10;
end

%% Save min solution (only needs to be run once) 
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
            vals1 = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS,:); % session 1 
            vals2 = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution(iS+1,:); % session 2
            vals12 = cat(1,vals1,vals2);
            meanVals12 = mean(vals12,1,'omitnan');
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolutionSub(subCount,:) = meanVals12;
        end
    end
end

% save mdlFit, fittedP, stats
matName = sprintf('ITPCfit_1_freeFreq_%s',dateStr); 
fullMatName = sprintf('%s/%s',figDir,matName); 
save(fullMatName,'mdlFit','-v7.3') % 'fittedP','stats'



