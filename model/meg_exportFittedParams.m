function meg_exportFittedParams
% exports fitted parameters as a CSV 

%% Settings
exportCSV = 1; 
csvDir = '/Users/kantian/Dropbox/Data/TANoise/TANoise-stats/data';
if ~exist(csvDir, 'dir')
    mkdir(csvDir)
end

mdlFitType = 'separate';
switch mdlFitType
    case 'separate' % fixed frequency
        % --- Precue T1 ---
        fvals(:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2)
        solutions(:,:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2)
        % --- Precue T2 ---
        fvals(:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2)
        solutions(:,:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2)
end
nPerms = size(fvals,1); 
nSubs = size(fvals,2); 
nParams = size(solutions,3); 

%% Restructure data into 10 sub x 2 sessions x 2 precues
clear fittedP
for iF = 1:numel(fitTypes)
    for iV = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution,3) % params
        subCount = 0;
        for iS = 1:2:20 % sessions
            subCount = subCount+1;
            clear val valCT1 valCT2
            % extract only parameters with minimum fval 
            [minVal1,idx1] = min( fvals(:,iS,1) );
            [minVal2,idx2] = min( fvals(:,iS,2) );

            valCT1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(idx1,iS:iS+1,iV);
            valCT2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(idx2,iS:iS+1,iV);
            val = squeeze(cat(4,valCT1,valCT2));
            fittedP.(fitTypes{iF}).(paramNames{iV})(subCount,:,:) = val;
        end
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
        % clear val valCT1 valCT2
        % valCT1 = mean(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).yhat(iS:iS+1,baselineToi),2,'omitnan');
        % valCT2 = mean(mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).yhat(iS:iS+1,baselineToi),2,'omitnan');
        % 
        % val = squeeze(cat(4,valCT1,valCT2));
        % fittedP.(fitTypes{iF}).baseline_yhat(subCount,:,:) = val;
    end
end

%% Prepare for export 
% Long format data for R
for iF = 1:numel(fitTypes)
    clear V
    fields = fieldnames(fittedP.(fitTypes{iF})); 
    switch fitTypes{iF}
        case 'linear2Hz'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope',...
                'amplitude','phase','freq'};
        case 'linear'
            tableHeaders = {'subject','session','precue',...
                'intercept','slope'};
    end
    count = 1;
    for iS = 1:10 % subjects
        for iC = 1:numel(cueLevel)
            for iSession = 1:2
                for iV = 1:numel(fields)
                    V.subject(count) = iS;
                    V.session(count) = iSession;
                    V.precue(count)  = iC;
                    V.(tableHeaders{iV+3})(count) = mdlFit.(fitTypes{iF}).(cueLevel{iC}).(fitLevel).solution(iP,iV,iS);
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

        csvName = sprintf('TANoise_ITPCFit_%s_%s',fitTypes{iF},dateStr);
        csvPath = sprintf('%s/%s.csv', csvDir, csvName);
        writetable(T,csvPath)
    end
end
