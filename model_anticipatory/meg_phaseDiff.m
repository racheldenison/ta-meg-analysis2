function meg_phaseDiff

% first run meg_fitAnticipatory1

%% Setup 
addpath(genpath('/Users/kantian/Dropbox/Software/watsons_u2'))

%% 
clear valCT1 valCT2
% test = fittedP.linear2Hz.phaserad(:,:,1); % all subjects and sessions, precue T1 
% testMean = circ_mean(test,[],2); 

%% Average sessions -> subjects  
valCT1 = circ_mean(fittedP.linear2Hz.phaserad(:,:,1),[],2); % average sessions 
valCT2 = circ_mean(fittedP.linear2Hz.phaserad(:,:,2),[],2); % average sessions 
phaseDiff = circ_dist(valCT1,valCT2); 

[stats.h stats.mu] = circ_mtest(phaseDiff, 0.1); 
% Throwing error: requirements for confidence levels not met 

%% Try sessions
clear fittedPS
count = 1;
for iS = 1:10
    for iSession = 1:2
        clear val
        val = fittedP.linear2Hz.phaserad(iS,iSession,:);

        % --- Format into sessions (20) x precue (2)
        fittedPS.linear2Hz.phaserad(count,:) = val;
        count = count + 1;
    end
end

clear valCT1 valCT2 phaseDiff stats
valCT1 = fittedPS.linear2Hz.phaserad(:,1); 
valCT2 = fittedPS.linear2Hz.phaserad(:,2); 
phaseDiff = circ_dist(valCT1,valCT2); 

% Plot to check 
figure 
for iS = 1:size(phaseDiff,1)
    rho = 1; sSize = 100; 
    % --- Differences --- 
    polarscatter(phaseDiff(iS),rho+4,sSize,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    hold on
    % --- Precue (T1) --- 
    polarscatter(valCT1(iS),rho+1,sSize,'filled','MarkerFaceColor',p.cueColors(1,:))   
    % --- Precue (T2) --- 
    polarscatter(valCT2(iS),rho+2,sSize,'filled','MarkerFaceColor',p.cueColors(2,:))  
    % --- Subject lines --- 
    polarplot([valCT1(iS) valCT2(iS)],[rho+1 rho+2],'color',[0.8 0.8 0.8],'LineWidth',0.7)
end
polarplot([circ_mean(phaseDiff) circ_mean(phaseDiff)],[0 rho+4],'k','LineWidth',2)
% 95CI of mean diff 
[stats.mu, stats.ul, stats.ll] = circ_mean(phaseDiff); 
polarplot([stats.ul stats.ul],[0 rho+4],'k','LineWidth',1)
polarplot([stats.ll stats.ll],[0 rho+4],'k','LineWidth',1)

polarplot([circ_mean(valCT1) circ_mean(valCT1)],[0 rho+1],'color',p.cueColors(1,:),'LineWidth',1)
polarplot([circ_mean(valCT2) circ_mean(valCT2)],[0 rho+2],'color',p.cueColors(2,:),'LineWidth',1)

%% Permutation test
nPerm = 1000; 
clear valCT1_Perm valCT2_Perm phaseDiff_Perm phaseDiffMean_Perm
for iP = 1:nPerm
    idx = randperm(size(phaseDiff,1)); 
    isOdd = mod(idx,2);
    for iS = 1:size(phaseDiff,1)
        if isOdd(iS) % if odd, swap precue labels 
            valCT1_Perm(iS,iP) = valCT2(iS); 
            valCT2_Perm(iS,iP) = valCT1(iS); 
        else % if even, 
            valCT1_Perm(iS,iP) = valCT1(iS); 
            valCT2_Perm(iS,iP) = valCT2(iS); 
        end
    end
    phaseDiff_Perm(:,iP) = circ_dist(valCT1_Perm(:,iP),valCT2_Perm(:,iP));
    phaseDiffMean_Perm(iP) = circ_mean(phaseDiff_Perm(:,iP)); 
end

% Plot histogram 
figure
set(gcf,'Position',[100 100 600 500])
subplot 211 
hold on 
meg_figureStyle
bins = -pi/2:0.1:pi/2; 
histogram(phaseDiffMean_Perm,bins)
xlabel('Phase difference (rad)')
ylabel('Count')
percentile = invprctile(phaseDiffMean_Perm,stats.mu);
txt = sprintf('%0.2f%',percentile); 
xline(stats.mu,'r',txt)

% histogram of absolute value of the diff? 
subplot 212 
hold on 
meg_figureStyle
bins = -pi/2:0.1:pi/2; 
histogram(abs(phaseDiffMean_Perm),bins)
xlabel('Abs phase difference (rad)')
ylabel('Count')
percentile = invprctile(abs(phaseDiffMean_Perm),abs(stats.mu));
txt = sprintf('%0.2f%',100-percentile); 
xline(abs(stats.mu),'r',txt)

%% Try Watson's U test (test whether two distributions are the identical)
[stats.p,stats.U2_H0] = watsons_U2_approx_p(valCT1,valCT2);
[stats.p, stats.U2_obs, stats.U2_H0] = watsons_U2_perm_test(valCT1,valCT2,10); 

%% 
% Try Raleigh's Z test 
phaseDiff = circ_dist(valCT1,valCT2); 
% Plot to check 
figure
for iS = 1:20
    theta = phaseDiff(iS); 
    rho = 1; sSize = 50; 
    polarscatter(theta,rho,sSize,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    hold on
end
[stats.pval, stats.z] = circ_rtest(valR_diff_min);

[stats.mu, stats.ul, stats.ll] = circ_mean(valR_diff_min'); 


% OLD 
% [stats.h,stats.p,stats.ci,stats.stats] = ttest(phaseDiff); 
% 
% meanAbsPhaseDiff = mean(abs(phaseDiff)); 
% 
% meanAbsPhaseDiff_deg = rad2deg(meanAbsPhaseDiff); 
% [stats.h,stats.p,stats.ci,stats.stats] = ttest(rad2deg(abs(phaseDiff))); 
% 
% [stats.h stats.mu] = circ_mtest(phaseDiff, 0); 