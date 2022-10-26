%% 
figure
subplot 311
title('RD, data ssvef workspace')
hold on 
t = -1500:5700; 
val = trigData; 
val = val(:,1,:); 
val = mean(val,3,'omitnan'); 
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlims = [-1000,5000]; 
xlim(xlims)
ylabel('fT')

subplot 312
title('KT, ebi mat')
hold on 
val = data; 
val = val(:,1,:); 
val = mean(val,3,'omitnan'); 
t = -2000:6000; 
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlims = [-1000,5000]; 
xlim(xlims)
ylabel('T')

subplot 313
title('KT, ebi mat, unit conversion')
hold on 
val = data; 
val = val(:,1,:); 
val = mean(val,3,'omitnan'); 
femtoUnit = 10e-15; 
val = val/femtoUnit; 
t = -2000:6000; 
plot(t,squeeze(val),'Color',[0.5 0.5 0.5],'LineWidth',2)
meg_figureStyle
xlims = [-1000,5000]; 
xlim(xlims)
ylabel('fT')