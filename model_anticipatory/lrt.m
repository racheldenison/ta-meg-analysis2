
% likelihood ratio test 

%% Load data and mdlFit 
% 2 Hz separate fits 
% load('ModelFit_SeparatePrecueT1T2_231006.mat')

% 2 Hz yoked fits 
load('ModelFit_YokedFixedFreq2Hz_PrecueT1T2_231002.mat')

%% Extract fvals (RMSE)
clear fvals 
for iF = 1:numel(fitTypes)
    for iS = 1:20
        clear fvalsPerm
        fvalsPerm = mdlFit.(fitTypes{iF}).session.fval(:,iS);
        [minVal,idx] = min( fvalsPerm );
        fvals.(fitTypes{iF}).sessions(iS) = fvalsPerm(idx);
    end
end

%% Average to subjects 
for iF = 1:numel(fitTypes)
    subCount = 0;
    for iS = 1:2:20 % sessions
        subCount = subCount+1;
        clear val
        val = fvals.(fitTypes{iF}).sessions(iS:iS+1);
        fvals.(fitTypes{iF}).sessions2(subCount,:) = val; % subjects(10) x sessions (2)
    end
    fvals.(fitTypes{iF}).subjects = mean(fvals.(fitTypes{iF}).sessions2,2,'omitnan'); 
end

%% Get degrees of freedom and likelihood ratio test 
% dof = number of restrictions in from full to restricted model 
dof = size(mdlFit.linear2Hz.session.solution,3) - size(mdlFit.linear.session.solution,3); 

for iS = 1:10
    uLogL = -fvals.linear2Hz.subjects(iS); % unrestricted model 
    rLogL = -fvals.linear.subjects(iS); % restricted model 
    [fvals.h(iS),fvals.pVal(iS),fvals.stat,fvals.cVal] = lratiotest(uLogL,rLogL,dof); 
end


