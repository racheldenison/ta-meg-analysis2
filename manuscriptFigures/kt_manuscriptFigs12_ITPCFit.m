% kt_manuscriptFigs12_ITPCFit

% August 2023
% Karen 

%% settings
saveFigs = 0; 
user = 'kantian'; % kantian karen 

%% Paths
addpath(genpath(pwd))

% Add circular stats toolbox 
% Berens (2009) https://www.jstatsoft.org/article/view/v031i10
addpath(sprintf('/Users/%s/Dropbox/Software/CircStat2012a',user)) 

%% load ITPC data 
load(sprintf('/Users/%s/Dropbox/github/ta-meg-analysis2/unused/groupA_ITPCspectrogram_byAtt.mat',user))
p = meg_params('TANoise_ITPCsession8');

%% Try with real data 
% y = A * sin(B(x + C)) + D; 
% amplitude = A 
% period = 2*pi/B 
% frequency = 1/period 
% From FFt we know the frequency to fit is 2 Hz , so B = 2/2*pi
% vertical shift = D 
% horizontal shift = C 

% --- Load ITPC data from above ---
sampling = 1:10:7001;
freq = 20; % frequency of interest 

% --- Fit settings ---
paddingBefore = 80; % ms before T1 
toi = abs(p.tstart)+p.eventTimes(1):abs(p.tstart)+p.eventTimes(2); % preCue:T1
toi = toi(1):toi(end)-paddingBefore; % select timing for fit 

% --- All trials --- 
% select ITPC at 20Hz > average to subjects > do fit 
sessionVals_all = squeeze(A.all.session(freq,:,:)); 
subVals_all = squeeze(A.all.subject(freq,:,:));  % frequency x time x subjects (from avg sessions)
groupVals_all = squeeze(mean(subVals_all,2)); % average across subjects 
% --- Precue T1 trials --- 
sessionVals_precueT1 = squeeze(A.cueT1.session(freq,:,:)); 
subVals_precueT1 = squeeze(A.cueT1.subject(freq,:,:)); 
groupVals_precueT1 = squeeze(mean(subVals_precueT1,2)); 
% --- Precue T2 trials --- 
sessionVals_precueT2 = squeeze(A.cueT2.session(freq,:,:)); 
subVals_precueT2 = squeeze(A.cueT2.subject(freq,:,:)); 
groupVals_precueT2 = squeeze(mean(subVals_precueT2,2)); 

% --- Prepare data as dataset for model ---
ds = dataset(toi', groupVals_all(toi)); 
vn = {'t','ITPC'}; 
ds = set(ds,'VarNames',vn);

%% --- Define model ---
% -- Change inputs here --
fitType = 'sineLinear'; % sineLinear, linear, sine 
nVars = 6; % number of free parameters in model 
% -------------------------

switch fitType
    case 'linear'
        switch nVars
            case 2
                modelfun = 'ITPC ~ (b1*t) + b2';
                % b1 is slope
                % b2 is intercept
                nVars = 2;
                beta0 = randn(nVars,1); % initial search
            otherwise
                error('fit not defined for this nVars')
        end
    case 'sine' % sine only, no linear component
        switch nVars
            % For 2 Hz frequency, at 1000 Hz sampling rate, then period = 500 samples
            % B = 2*pi/period
            % B = 2*pi/(500) = 0.0126
            case 2 
                % b1 phase shift
                % b2 vertical shift
                % fixed: amplitude and frequency 
                modelfun = 'ITPC ~ 0.1 * sin( 0.0126 * (t + b1)) + b2';
                beta0 = randn(nVars,1);
            case 3
                % b1 = amplitude
                % b2 = phase shift
                % b3 = vertical shift
                % fixed: frequency 
                modelfun = 'ITPC ~ b1 * sin( 0.0126 * (t + b2)) + b3';
                beta0 = randn(nVars,1);
                beta0(1) = 0.1;
                beta0(2) = 255;
                beta0(3) = 0.34;
            case 4
                % b1 = amplitude
                % b2 = phase shift
                % b3 = vertical shift
                % b4 = freq*pi/(500)
                % fixed: nothing 
                modelfun = 'ITPC ~ b1 * sin( b4 * (t + b2)) + b3';
                beta0 = randn(nVars,1);
                beta0(1) = 0.05;
                beta0(2) = 255;
                beta0(3) = 0.34;
                beta0(4) = 0.0126;
        end
    case 'sineLinear'
        switch nVars
            case 4
                % b1 = phase shift
                % b2 = vertical shift
                % b3 = slope
                % b4 = freq*pi/(500)
                modelfun = 'ITPC ~ (b3*t) + ( 0.1 * sin( 0.0126 * (t + b1)) + b2) + b4';
                beta0 = randn(nVars,1);
            case 6 
                % b1 = amplitude
                % b2 = phase shift
                % b3 = vertical shift
                % b4 = freq*pi/(500)
                paramVarNames = {'b1','b2','b3','b4','b5','b6'};
                paramNames = {'amplitude','phase','vertical','freqB','slope','intercept'};
                fixedParams = {};
                modelfun = 'ITPC ~ (b5*t) + ( b1 * sin( b4 * (t + b2)) + b3) + b6';
                beta0 = randn(nVars,1);
                beta0(1) = 0.5;
                beta0(2) = 255; % phase in deg 
                beta0(3) = 0.34;
                beta0(4) = 0.0126;
                beta0(5) = 1e-06; % 1.1e-06
                beta0(6) = 0.5; 
        end
end

% --- Fit model to data (group, for visualization) ---
clear mdl_all_group mdl_precueT1_group mdl_precueT2_group
mdl_all_group = fitnlm(toi',groupVals_all(toi),modelfun,beta0); % returns estimated coefficients 
mdl_precueT1_group = fitnlm(toi',groupVals_precueT1(toi),modelfun,beta0);
mdl_precueT2_group = fitnlm(toi',groupVals_precueT2(toi),modelfun,beta0);

% --- Do fit (on subjects) ---
clear estCoeffs mdl_all mdl_precueT1 mdl_precueT2
for i = 1:size(subVals_all,2)
    mdl_all = fitnlm(toi',subVals_all(toi,i),modelfun,beta0); % returns estimated coefficients
    mdl_precueT1 = fitnlm(toi',subVals_precueT1(toi,i),modelfun,beta0);
    mdl_precueT2 = fitnlm(toi',subVals_precueT2(toi,i),modelfun,beta0);
    % save to matrix
    estCoeffs.all(:,i) = mdl_all.Coefficients.Estimate;
    estCoeffs.precueT1(:,i) = mdl_precueT1.Coefficients.Estimate;
    estCoeffs.precueT2(:,i) = mdl_precueT2.Coefficients.Estimate;
end

% --- ttests on fitted model params (Precue T1 vs Precue T2) --- 
clear F
for i = 1:nVars
    F(i).paramNames = paramNames{i};
    switch paramNames{i}
        case {'phase'} % circular stats 
            clear val1 val2
            val1 = circ_ang2rad(estCoeffs.precueT1(i,:)); 
            val2 = circ_ang2rad(estCoeffs.precueT2(i,:)); 
            F(i).p = circ_wwtest(val1,val2);
            % --- polar plot --- 
            figure
            set(gcf,'Position',[100 100 400 400])
            subplot 211 
            hold on
            circ_plot(mean(val1),'pretty','bo',true,'linewidth',2,'color',p.cueColors(1,:))
            subplot 212
            hold on
            circ_plot(mean(val2),'pretty','bo',true,'linewidth',2,'color',p.cueColors(2,:))
            % -- plot error bars --- 
        case {'amplitude','vertical','freqB','slope','intercept'}
            [F(i).h, F(i).p, F(i).ci, F(i).stats] = ttest( estCoeffs.precueT1(i,:) , estCoeffs.precueT2(i,:) );
        otherwise 
            error('param for ttest unrecognized') 
    end
end

%% --- Figure (subjects) ---
figure
set(gcf,'Position',[100 100 800 800])

% fitType = 'linear';
for i = 1:size(subVals_all,2)
    subplot (size(subVals_all,2)/2,2,i)
    hold on

    % --- Plot data ---
    plot(p.t, subVals_precueT1(:,i),'LineWidth',3,'Color',p.cueColors(1,:))
    plot(p.t, subVals_precueT2(:,i),'LineWidth',3,'Color',p.cueColors(2,:))

    % --- Plot fit ---
    switch fitType
        case 'linear'
            % --- Precue T1 fit ---
            yhat = (estCoeffs.precueT1(1,i)*toi) + estCoeffs.precueT1(2,i);
            plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
            % --- Precue T2 fit ---
            yhat = (estCoeffs.precueT2(1,i)*toi) + estCoeffs.precueT2(2,i);
            plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
        case 'sine'
            switch nVars
                case 2
                    % --- Precue T1 fit ---
                    yhat = 0.1 * sin( 0.0126 * (toi + estCoeffs.precueT1(1,i))) + estCoeffs.precueT1(2,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
                    % --- Precue T2 fit ---
                    yhat = 0.1 * sin( 0.0126 * (toi + estCoeffs.precueT2(1,i))) + estCoeffs.precueT2(2,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
                case 3
                    % --- Precue T1 fit ---
                    yhat = estCoeffs.precueT1(1,i) * sin( 0.0126 * (toi + estCoeffs.precueT1(2,i))) + estCoeffs.precueT1(3,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
                    % --- Precue T2 fit ---
                    yhat = estCoeffs.precueT2(1,i) * sin( 0.0126 * (toi + estCoeffs.precueT2(1,i))) + estCoeffs.precueT2(2,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
            end
        case 'sineLinear'
            switch nVars
                case 2
                    % --- Precue T1 fit ---
                    yhat = (estCoeffs.precueT1(3,i)  * toi) + (0.1 * sin( 0.0126 * (toi + estCoeffs.precueT1(1,i))) + estCoeffs.precueT1(2,i)) + estCoeffs.precueT1(4,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
                    % --- Precue T2 fit ---
                    yhat = (estCoeffs.precueT2(3,i)  * toi) + (0.1 * sin( 0.0126 * (toi + estCoeffs.precueT2(1,i))) + estCoeffs.precueT2(2,i)) + estCoeffs.precueT2(4,i);
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
                case 6
                    t = p.t(toi); 
                    % --- Precue T1 fit ---
                    bFit = estCoeffs.precueT1; 
                    yhat = ( bFit(5,i) *t) + ( bFit(1,i) * sin( bFit(4,i) * (t + bFit(2,i))) + bFit(3,i)) + bFit(6,i); 
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
                    % --- Precue T2 fit ---
                    bFit = estCoeffs.precueT2; 
                    yhat = ( bFit(5,i) *t) + ( bFit(1,i) * sin( bFit(4,i) * (t + bFit(2,i))) + bFit(3,i)) + bFit(6,i); 
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
                case 'groupFit'
                    % --- Precue T1 fit ---
                    yhat = (4.6796e-05  * toi) + (0.1 * sin( 0.0126 * (toi + estCoeffs.precueT1(1,i))) + estCoeffs.precueT1(2,i)) + 0.20348;
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(1,:))
                    % --- Precue T2 fit ---
                    yhat = (4.6796e-05  * toi) + (0.1 * sin( 0.0126 * (toi + estCoeffs.precueT2(1,i))) + estCoeffs.precueT2(2,i)) + 0.20348;
                    plot(p.t(toi),yhat,'LineWidth',1,'Color',p.cueColors(2,:))
            end

    end

    % --- Plot event lines ---
    for i = 1:numel(p.eventTimes)
        xline(p.eventTimes(i),'Color',[0.5 0.5 0.5],'LineWidth',1)
    end

    % --- Format --- 
    meg_figureStyle
    xlim([-100 2400])
end

% --- Plot fit --- 
% p1_fit = plot(p.t(toi),y1_est,'Color',[0.3 0.3 0.3],'LineWidth',2); 

%% --- Figure (group) ---
figure
set(gcf,'Position',[100 100 600 400])
hold on 

% --- Plot fit --- 
p1_fit = plot(p.t(toi),y1_est,'Color',[0.3 0.3 0.3],'LineWidth',2); 

% --- Downsample ---
groupVals = groupVals(sampling); 





















