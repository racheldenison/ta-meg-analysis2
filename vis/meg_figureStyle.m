function meg_figureStyle()

% adjusts fig axis and text styling
box off
set(gca,'TickDir','out');
ax = gca;
ax.LineWidth = 1.5;
ax.XColor = 'black';
ax.YColor = 'black';

smlFont = 14;
bigFont = 24; 

ax.FontSize = bigFont;
ax.FontName = 'Helvetica-Light'; 
 
ax.XAxis.FontSize = smlFont;
ax.YAxis.FontSize = smlFont;

ax.LabelFontSizeMultiplier = bigFont/smlFont; 








