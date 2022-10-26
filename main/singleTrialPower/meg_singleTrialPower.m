% Single trial analysis 
% September 15, 2022

for i = 1:20 
    val = squeeze(groupA(i).all.tfPowsTrials); 
    val = mean(val,3,'omitnan'); 
    groupA(i).all.singleTrialPow(:,:,i) = val; 
end


%% 
val = []; 
t = p.tstart:p.tstop; 
toi = -100:2000; 
toi = -500:5000;
trials = 1:3; 
% trials = 50:100; 

tOIIdx = find(t==toi(1)):find(t==toi(end)); 

figure
subplot 611
groundTruth = zeros(size(toi));
period = (1/20)*1000; % ms 
for i = 1:period:size(toi,2)-period
    groundTruth(i:i+period/2) = 1; 
end
plot(toi,groundTruth)
meg_figureStyle
xlim([toi(1) toi(end)])

for iCh = 1:5 
subplot (6,1,iCh+1)
i = 1; 
val = groupA(i).all.tfPowsTrials(trials,:,tOIIdx,iCh); % trial x time x session
val = squeeze(val); 
% imagesc(val) 
plot(toi,val)
colormap(gray)
meg_figureStyle
xlim([toi(1) toi(end)])
end

