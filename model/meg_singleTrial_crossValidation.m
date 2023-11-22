function meg_singleTrial_crossValidation(expt, sessionDir, user)
% function meg_singleTrial_crossValidation(expt, sessionDir, user)
% split half cross validation per subject 

analStart = tic; % for timing check 

%% Setup 
% sessionDir = 'R0817_20171212'; 
% expt = 'TANoise'; 

cluster = 0; % using cluster? 0=local, 1=cluster 

[sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt); 

sessionIdx = find(strcmp(sessionNames,sessionDir)); 

condNames = {'all'}; 
nChannels = 5; 

user = 'kantian'; 
workingDir = sprintf('/Users/%s/Dropbox/github/ta-meg-analysis-model/model_anticipatory',user); 
addpath(genpath(workingDir)) % ta-meg-analysis-model/model_anticipatory
addpath(genpath('../../ta-meg-analysis2')) 

figFormat = 'png';

dateStr = datetime('now','TimeZone','local','Format','yyMMdd');
figDir = sprintf('%s/singleTrialITPC/%s',workingDir,dateStr); 
if ~exist(figDir, 'dir')
    mkdir(figDir)
end

% --- Add circular stats toolbox ---
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

%% Analysis settings 
% --- Cross validation split half setttings ---
nPermCV = 5; % 10, 100 
% --- Model fit start coefficient perms ---
nPerms = 3; % 10, 100

% --- MEG settings --- 
p = meg_params('TANoise_ITPCsession8');

% --- Timing and ITPC settings --- 
tT = p.tstart:p.tstop;

taper          = 'hanning';
foi            = 20; % 100;
t_ftimwin      = 10 ./ foi;
toiT            = p.tstart/1000:0.01:p.tstop/1000; % toi for (T)rial 
tfAmps = []; tfAmpsAvg = [];
tfPows = []; tfPowsAvg = [];

padTotal = ceil(p.trialTime/p.fSample);
padPre = ceil((padTotal*p.fSample-p.trialTime)/2);
padPost = floor((padTotal*p.fSample-p.trialTime)/2);
toiPad = (p.tstart-padPre)/1000:0.01:(p.tstop+padPost)/1000;
tPad = p.tstart-padPre:p.tstop+padPost;
xtick = 1:80:numel(toiPad);
ytick = 10:10:numel(foi);
xlims = [size(toiT,1),size(toiT,2)];
Fsample = p.fSample;
width = 8;

% --- Model fitting settings 
freq = 2; % Hz, fixed frequency for model fit 
fitLevel = 'session'; 
% anticipatory toi 
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore;
tIdx = toi+1; % time index
t = p.t(tIdx)+1; % trial relative time 
Fs = 1000; 

baselineToi = abs(p.tstart)+p.eventTimes(1)-100:abs(p.tstart)+p.eventTimes(1); 

%% --- Fit settings --- 
paramNames = {'intercept','slope','amplitude','phase'};
nVars = numel(paramNames); 

% === SETTINGS ===
fitLevel = 'session'; % 'session' 'group'
cueLevel = {'all'}; % 'cueT1','cueT2' 'all'
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

mdlFit.lb = lb; 
mdlFit.ub = ub; 

mdlFit.hardlb = mdlFit.lb;
mdlFit.hardub = mdlFit.ub;

% --- Nonlinear constraints --- 
nonlcon = []; 

% --- Options --- 
options = optimoptions('fmincon','Display','iter'); 

%% Run analysis
% Load data 
% clear group A4
% 
% [A, D, selectedChannels, I, B] = meg_runAnalysis(expt, sessionNames{iS}, 'karen');
% group(iS).A = A;
% group(iS).D = D;
% group(iS).selectedChannels = selectedChannels;
% group(iS).B = B;

% edit to be dynamic path 
% groupA: MEG, 20 Hz phase angle: trials (384) x ch (5) x foi (1) x time (7001)
% groupB 
load(sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/data/singleTrialPower_allTrials_s20.mat',user)); 

%% Calculate ITPC from split training testing data
% A3 for analysis of split training testing data
clear A3
phaseAngle = squeeze(groupA(sessionIdx).all.phaseAngle); % trials (384) x ch (5) x time (7001) 
nTrials = size(phaseAngle,1);
nChannels = size(phaseAngle,2); 

for iPerm = 1:nPermCV % permute splits for cross validation
    clear ITPCtraining ITPCtesting
    idxRand = randperm(nTrials);

    nTrialsTest = nTrials/2; % for split half
    nTrialsTrain = nTrials - nTrialsTest;

    trainingTrials = idxRand(1:nTrialsTest);
    testingTrials = idxRand(nTrialsTest+1:nTrials);

    % Get phase angles by training & testing splits 
    for iCh = 1:nChannels
        for iT = 1:2 % training, then testing
            clear vals idxTrials
            if iT==1 % training
                idxTrials = trainingTrials; 
            elseif iT==2 % testing 
                idxTrials = testingTrials; 
            end
            vals = squeeze(phaseAngle(idxTrials,iCh,:)); % trials (192) x time (7001) 

            % Calculate ITPC from phase angles 
            itpc = squeeze(abs(mean(exp(1i*vals),1,'omitnan')));

            if iT==1
                ITPCtraining(:,iCh) = itpc; % f x t x  ch
            elseif iT==2
                ITPCtesting(:,iCh) = itpc;
            end
        end
    end
    
    % Save ITPC by split halves, for model fit 
    A3.all.ITPCMean.training(:,iPerm) = nanmean(ITPCtraining,2); % average channels
    A3.all.ITPCMean.testing(:,iPerm) = nanmean(ITPCtesting,2); 
    % Save idx of training and testing trials 
    A3.all.trialsIdx.training(:,iPerm) = trainingTrials; 
    A3.all.trialsIdx.testing(:,iPerm) = testingTrials; 

    %% --- Define initial search array ---
    nGrain = 100;
    clear x0s
    x0s = NaN([nVars,nGrain]);
    paramNames = {'intercept','slope','amplitude','phase'};
    for iStart = 1:nVars
        x0s(iStart ,:) = linspace(mdlFit.lb(iStart),mdlFit.ub(iStart),nGrain); % initial coefficients for search
    end

    %% Both model fits (linear, linear+2Hz)
    clear x0 x0Perm
    for iC = 1:numel(cueLevel)
        clear dataFit
        dataFit = A3.(cueLevel{iC}).ITPCMean.training(tIdx,iPerm)'; % 1 x time (971), fit data on training partition
        dataTest = A3.(cueLevel{iC}).ITPCMean.testing(tIdx,iPerm)';

        for iF = 1:numel(fitTypes)
            for iP = 1:nPerms % iP is optimization permutation start grid
                if ~cluster
                    txt = sprintf('Fitting %s, permutation %d of %d...',sessionDir,iP,nPerms);
                    disp(txt)
                end
                % --- Randomly pick starting coefficients from search grid ---
                for iV = 1:nVars
                    vals = randperm(nGrain);
                    idx(iV) = vals(1);
                    x0(iV) = x0s(iV,idx(iV));
                end

                switch fitTypes{iF}
                    case 'linear2Hz'
                        if ~cluster
                            disp('Fitting linear + 2Hz ...')
                        end
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
                        if ~cluster
                            disp('Fitting linear ...')
                        end
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
                % mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).paramNames = paramNames;
                % -- Save fixed freq ---
                % mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fixedFreq = freq; % hz
                % --- Save fitted parameters ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,:) = solution;
                % --- Save error ---
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(iP) = fval;
                % --- Calculate r squared --- 
                clear E yhat rsq
                [E,yhat] = meg_objectiveFunction1(solution,dataFit,t,Fs,paramNames,fitTypes{iF});
                rsq = calculateRSQ(dataTest,yhat,1);
                mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).rsq(iP) = rsq; 
                % --- Save predicted y ---
                % mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).yhat(iP,iS,:) = yhat;
            end % end the search permutation 

            % --- Find permutation with minimum fval --- 
            [minVal,idx] = min(  mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).fval(:)  );
            fittedX = squeeze(mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(idx,:));
            % --- Save min fitted params per permutation of starting coefficients
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minSolution = fittedX;
            % --- Calculate R^2 ---
            %of fitted model on the test data
            switch fitTypes{iF}
                case 'linear2Hz' 
                    paramNames = {'intercept','slope','amplitude','phase'};
                case 'linear'
                    paramNames = {'intercept','slope'};
            end
            [E,yhat] = meg_objectiveFunction1(fittedX,dataFit,t,Fs,paramNames,fitTypes{iF},freq);
            rsq = calculateRSQ(dataFit,yhat,1);
            % txt = sprintf('Rsq = %0.2f...',rsq);
            % disp(txt)
            % --- Save R^2 ---
            mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).minRsq = rsq;
            
            %% Check fit on test partition --> R^2
            clear rsq
            rsq = calculateRSQ(dataTest,yhat,1);

            %% Plot for testing (only first 10 perms) 
            if iPerm<=10
                if iPerm==1 && iF==1
                    figure
                    set(gcf,'Position',[100 100 500 900])
                    sgtitle(und2space(sessionDir))
                end
                subplot (5,2,iPerm)
                hold on
                meg_figureStyle
                plot(dataTest)
                plot(yhat)
                plot(dataFit)
                xlabel('Time (ms)')
                ylabel('ITPC')
                if iPerm==10 && iF==2
                    figTitle = sprintf('%s_TANoise_CrossValidation_first10Perms',sessionDir);
                    saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))
                end
            end

            %% Save all
            A4.(fitTypes{iF}).(cueLevel{iC}).rsq(iPerm) = rsq;
        end
    end
end

%% Plot histogram
figure
set(gcf,'Position',[100 100 400 370])
sgtitle(und2space(sessionDir))

edges = -1:0.1:1; 
subplot 121
hold on
meg_figureStyle
histogram(A4.linear.all.rsq,edges)
histogram(A4.linear2Hz.all.rsq,edges)
xlabel('R squared')
ylabel('Count')
l = legend('Linear','Linear + 2Hz');
l.FontSize = 9;

subplot 122
hold on
meg_figureStyle
bar(1,mean(A4.linear.all.rsq),'FaceColor',[0 0.4470 0.7410])
errorbar(1,mean(A4.linear.all.rsq), std(A4.linear.all.rsq)/sqrt(20),'CapSize',0,'Color','k','LineWidth',1.5)
bar(2,mean(A4.linear2Hz.all.rsq),'FaceColor',[0.8500 0.3250 0.0980])
errorbar(2,mean(A4.linear2Hz.all.rsq), std(A4.linear2Hz.all.rsq)/sqrt(20),'CapSize',0,'Color','k','LineWidth',1.5)
xlim([0 3])
xticks([1 2])
xticklabels({'Linear','Linear + 2Hz'})
ylabel('R squared')

% --- Save fig ---
figTitle = sprintf('%s_TANoise_CrossValidation',sessionDir);
saveas(gcf,sprintf('%s/%s.%s', figDir, figTitle, figFormat))

A4.analTiming = toc(analStart); % elapsed timing (s) 

%% Save A4
filename = sprintf('%s/CV_%s_%s',figDir,sessionDir,dateStr);
save(filename,'A4','mdlFit','A3','-v7.3')




