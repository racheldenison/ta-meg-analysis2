function megAvgFittedParams(mdlFit,paramNames,mdlFitType)
% function megAvgFittedParams
% averages fitted parameters to subject level 

%% Setup 
switch mdlFitType
    case 'separate' % fixed frequency
        idxPhase = [find(contains(paramNames,'phase'))];
        % --- Precue T1 ---
        fvals(:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2) 
        solutions(:,:,:,1) = mdlFit.linear2Hz.cueT1.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2) 
        % --- Precue T2 ---
        fvals(:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).fval; % permutations (100) x subjects (20) x precue (2) 
        solutions(:,:,:,2) = mdlFit.linear2Hz.cueT2.(fitLevel).solution; % permutations (100) x subjects (20) x fitted parameters (4) x precue (2) 
    case 'yoked_fixedFreq'
        idxPhase = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))];
        fvals = mdlFit.linear2Hz.(fitLevel).fval;
        solutions = mdlFit.linear2Hz.(fitLevel).solution; % permutations (100) x subjects (20) x parameters (8) 
    case 'yoked'
        idxPhase = [find(contains(paramNames,'phase1')) find(contains(paramNames,'phase2'))];
        idxFreq = [find(contains(paramNames,'freq'))]; 
        fvals = mdlFit.linear2Hz.(fitLevel).fval;
        solutions = mdlFit.linear2Hz.(fitLevel).solution; % permutations (100) x subjects (20) x parameters (9) 
end
nPerms = size(fvals,1); 
nSubs = size(fvals,2); 
nParams = size(solutions,3); 
cueLevel = {'cueT1','cueT2'};

%% --- Extract fitted phases  --- 
for iS = 1:nSubs
    figure
    set(gcf,'Position',[100 100 600 200])
    switch mdlFitType
        case 'separate'
            freq = 2; % Hz 
            % --- Precue T1 ---
            [minVal1,idx1] = min( fvals(:,iS,1) );
            valR1 = solutions(:,iS,idxPhase,1); % radian space 
            valR1_min = valR1(idx1);
            % --- Precue T2 ---
            [minVal2,idx2] = min( fvals(:,iS,2) );
            valR2 = solutions(:,iS,idxPhase,2); % radian space 
            valR2_min = valR2(idx2); 
            % --- Phase difference ---
            valR_diff = abs(circ_dist(valR1,valR2));
            valR_diffMin = abs(circ_dist(valR1_min,valR2_min)); 
        case 'yoked_fixedFreq'
            freq = 2; % Hz 
            [minVal,idx] = min( fvals(:,iS) );
             % --- Precue T1 ---
            valR1 = solutions(:,iS,idxPhase(1)); % radian space 
            valR1_min = valR1(idx);
            % --- Precue T2 ---
            valR2 = solutions(:,iS,idxPhase(2)); % radian space 
            valR2_min = valR2(idx); 
        case 'yoked'
            [minVal,idx] = min( fvals(:,iS) );
            freq = solutions(idx,iS,idxFreq); % Hz
            % --- Precue T1 ---
            valR1 = solutions(:,iS,idxPhase(1)); % radian space
            valR1_min = valR1(idx);
            % --- Precue T2 ---
            valR2 = solutions(:,iS,idxPhase(2)); % radian space
            valR2_min = valR2(idx);
    end
    % --- Phase difference ---
    valR_diff = abs(circ_dist(valR1,valR2));
    valR_diffMin = abs(circ_dist(valR1_min,valR2_min));
end

%% 
clear fittedP
iP = 1;
for iF = 1:numel(fitTypes)
    for iV = 1:size(mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution,2) % params
        subCount = 0;
        for iS = 1:2:20 % sessions
            subCount = subCount+1;
            clear val valCT1 valCT2
            valCT1 = mdlFit.(fitTypes{iF}).(cueLevel{1}).(fitLevel).solution(iP,iV,iS:iS+1);
            valCT2 = mdlFit.(fitTypes{iF}).(cueLevel{2}).(fitLevel).solution(iP,iV,iS:iS+1);
            val = squeeze(cat(4,valCT1,valCT2));
            fittedP.(fitTypes{iF}).(paramNames{iV})(subCount,:,:) = val;
        end
    end
end

