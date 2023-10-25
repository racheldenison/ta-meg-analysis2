% plots power law cartoon
f = 1:100; 
p = 1./f; 
figure
set(gcf,'Position',[100 100 500 200])

subplot 121
hold on 
meg_figureStyle
grey = 0.7*[1 1 1]; 
width = 5; 
plot(f,p,'LineWidth',width,'Color',grey)
ylabel('Amplitude')
xlabel('Frequency')
ax = gca; 
set(ax ,'Layer', 'Top')

subplot 122 
hold on 
meg_figureStyle
plot(log(f),log(p),'LineWidth',width,'Color',[grey])
ylabel('Log (amplitude)')
xlabel('Log (frequency)')
ax = gca; 
set(ax ,'Layer', 'Top')

saveas(gcf,'powerLaw.jpg')



