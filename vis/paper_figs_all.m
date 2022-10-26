% paper_figs_all.m

% Paper figures for Nat Neuro submission

%% Will
% home = '~'
% map = load('~/Google Drive/MATLAB/utilities/MyColorMaps.mat');
% load('~/Google Drive/Will - Confidence/Analysis/optimizations/attention3_MASTER.mat');
% datadir = '/Users/will/Google Drive/Will - Confidence/Data/attention3_ab';

%% directories
home = '~';
% home = '/Local/Users/denison/';
datadir = sprintf('%s/Google Drive/Shared/Projects/Will - Confidence/Data/attention3_ab',home);
figDir = sprintf('%s/Google Drive/Shared/Projects/Attention_Uncertainty/Figures',home);

%% add paths
addpath(genpath(sprintf('%s/Google Drive/Shared/Code/utilities',home)));
addpath(sprintf('%s/Google Drive/NYU/Projects/Confidence_Distributions/Code/confidence/helper_functions',home))

%% load data and maps
load(sprintf('%s/Google Drive/Shared/Projects/Will - Confidence/Analysis/optimizations/attention3_MASTER.mat',home))
map = load('MyColorMaps.mat');
ss = load('fig2sorted_subject_numbers.mat');
reorder_subjects = ss.subNums(:,1);

%% MAIN TEXT PANELS
%% distributions
hDist = fill_criteria_plots;

hDist(3).Children(2).XLabel.String = 'Internal measurement of orientation (°)';
hDist(4).Children(1).String = {'Category 1','Category 2'};

%% tetrapus
hPus = static_model_diagrams;

for i = 1:numel(hPus.Children)
    title(hPus.Children(i),'')
end

%% accuracy, confidence, RT % can we not do bars plz % i <3 barz
ax.col = 'depvar';
ax.row = 'slice';
ax.fig = 'none';
tasks = {'B'};

% 'c' = validity, 'c_s' = validity x orientation
show_data_or_fits('root_datadir', datadir, 'nBins', 7, 'depvars', {'tf', 'g', 'rt'}, 'axis', ax, ...
    'slices', {'c', 'c_s'}, 'tasks', tasks, 'symmetrify', false, 'errorbarwidth', 1.6, 'linewidth', 1.5,...
    'meanlinewidth', 1.5, 'mean_color', 'k', 'gutter', [.05 .1], 'margins', [.1 .1 .1 .1], ...
    's_labels', [-10 -5 -2 0 2 5 10], 'label_s_bin_centers', true, 'show_legend', true, ...
    'panel_size', [200 250]);

h = gcf;
h.Children(7).YLim = [.5 .8];
% h.Children(4).YLim = [1 3.5];

ylabel(h.Children(1), 'Reaction time (s)');
ylabel(h.Children(2), 'Reaction time (s)');
ylabel(h.Children(3), 'Mean confidence');
ylabel(h.Children(4), 'Mean confidence');
ylabel(h.Children(6), 'Prop. correct');
ylabel(h.Children(7), 'Prop. correct');

xlabel(h.Children(1),'')
xlabel(h.Children(2),'')
xlabel(h.Children(3),'')
xlabel(h.Children(4),'')
xlabel(h.Children(6), 'Orientation (°)');
xlabel(h.Children(7), 'Cue validity');

h.Children(2).XTickLabel = {'V','N','I'};
h.Children(4).XTickLabel = {'V','N','I'};
h.Children(7).XTickLabel = {'V','N','I'};

hPerf = figure('Position', h.Position);
copyobj(h.Children, hPerf)

%% mean button press + fits and model comparison: Choice + confidence, Bayes and Fixed only
ax.col = 'model';
ax.row = 'slice';
ax.fig = 'none';
tasks = {'B'};

% Choice + confidence: Fixed, Bayes, Lin, Quad
ah = show_data_or_fits('models', modelmaster([2 7]), 'root_datadir', datadir, 'depvars', {'resp','tf','Chat', 'g'}, 'axis', ax, ...
    'slices', {'c_s'}, 'tasks', tasks, 'symmetrify', false, 'errorbarwidth', 2.3, 'linewidth', 1.5,...
    'meanlinewidth', 1.5, 'gutter', [.01 .08], 'margins', [.1 .1 .1 .1], ...
    's_labels', [-9 -5 -2 0 2 5 9], 'show_legend', true, 'nPlotSamples', 20, 'nFakeGroupDatasets', 1000,...
    'MCM', 'loopsis', 'MCM_size', .5, 'nRespSquares', 6, 'ref_model', 2, 'label_s_bin_centers', true, ...
    'sort_subjects', false, 'keep_subjects_sorted', true, 'show_subject_names', true, 'sort_model', 1, ...
    'reorder_subjects', reorder_subjects);

%     'sort_subjects', true, 'keep_subjects_sorted', true, 'show_subject_names', true, 'sort_model', 1);

h = gcf;
title(h.Children(5),'Bayesian')
ylabel(h.Children(5),'Mean response')
xlabel(h.Children(5),'Orientation (°)')
xlabel(h.Children(3),'')

hFits = figure('Position', h.Position);
copyobj(h.Children, hFits)

%% non-parametric choice boundaries
hBound = nonparametric_boundary_plots;

set(gca,'LineWidth',1)
xlabel('Orientation uncertainty (°)')
ylabel('Choice decision boundary (°)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN TEXT FIGURES
%% FIGURE 1
% fig1H = figure('Name','FIG 1','Position',[20 745 800 600]);
fig1H = figure('Name','FIG 1','Units','inches','Position',[1 10 11 8.5]);

% panel a: stimulus distributions
ax1H(1) = copyobj(hDist(4).Children(2), fig1H);

% panel b: trial sequence
% (leave blank)

% panel c: performance bar graphs
srcAx = [2 4 7];
destAx = [2 3 4];
for i = 1:numel(srcAx)
    ax1H(destAx(i)) = copyobj(hPerf.Children(srcAx(i)), fig1H);
end

% panel d: filled criteria plots
srcAx = 1:6;
destAx = 5:10;
for i = 1:numel(srcAx)
    ax1H(destAx(i)) = copyobj(hDist(3).Children(srcAx(i)), fig1H);
end

% panel e: tetrapus plots
srcAx = 1:2;
destAx = 11:12;
for i = 1:numel(srcAx)
    ax1H(destAx(i)) = copyobj(hPus.Children(srcAx(i)), fig1H);
end

% set axis properties
for i = 1:numel(ax1H)
    ax1H(i).Units = 'inches';
    ax1H(i).FontSize = 14;
    ax1H(i).LineWidth = 1;
end

% set positions of subplots
ax1H(1).Position = [.5 6.5 3 1.5];

ax1H(2).Position = [9.1 6.5 1.5 1.5];
ax1H(3).Position = [6.8 6.5 1.5 1.5];
ax1H(4).Position = [4.5 6.5 1.5 1.5];

ax1H(5).Position = [7 2.5 2.25 .9];
ax1H(6).Position = [4.5 2.5 2.25 .9];
ax1H(7).Position = [7 3.5 2.25 .9];
ax1H(8).Position = [4.5 3.5 2.25 .9];
ax1H(9).Position = [7 4.5 2.25 .9];
ax1H(10).Position = [4.5 4.5 2.25 .9];

ax1H(11).Position = [7 .8 2.25 .9];
ax1H(12).Position = [4.5 .8 2.25 .9];


%% FIGURE 2
% fig2H = figure('Name','FIG 2','Position',[40 725 800 600]);
fig2H = figure('Name','FIG 2','Units','inches','Position',[1 10 11 7.5]);

% panel a: performance as a function of orientation
srcAx = [1 3 6];
destAx = [1 2 3];
for i = 1:numel(srcAx)
    ax2H(destAx(i)) = copyobj(hPerf.Children(srcAx(i)), fig2H);
end

% panel b: mean response with fits
srcAx = [3 5];
destAx = [4 5];
for i = 1:numel(srcAx)
    ax2H(destAx(i)) = copyobj(hFits.Children(srcAx(i)), fig2H);
end

% panel c: model comparison bars
ax2H(6) = copyobj(hFits.Children(2), fig2H);
ax2H(6).Position = [.7 .2 .25 .2];

% set axis properties
for i = 1:numel(ax2H)
    ax2H(i).Units = 'inches';
    ax2H(i).FontSize = 14;
    ax2H(i).LineWidth = 1;
end

% set positions of subplots
ax2H(1).Position = [7.5 5 2.25 2];
ax2H(2).Position = [4.5 5 2.25 2];
ax2H(3).Position = [1.5 5 2.25 2];

ax2H(4).Position = [4.25 1.5 2.5 2.25];
ax2H(5).Position = [1.5 1.5 2.5 2.25];

ax2H(6).Position = [7.75 1.75 2.75 1.5];

%% FIGURE 3
fig3H = figure('Name','FIG 3','Position',[60 705 560 420]);
copyobj(hBound.Children, fig3H)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPLEMENT PANELS
%% choice boundaries illustration
hCBT = choice_boundary_plots_theoretical;

h = gcf;
xlabel(h.Children(2),'Orientation uncertainty (°)')
ylabel(h.Children(2),'Choice boundary (°)')

%% fitted choice boundaries for different models, one observer
hCBExample = plot_choicebounds_singlesubject;

h = gcf;
title(h.Children(2),'Observer 5')
xlabel(h.Children(2),'Orientation uncertainty (°)')
ylabel(h.Children(2),'Choice boundary (°)')

%% mean button press + fits and model comparison: Choice + confidence
ax.col = 'model';
ax.row = 'slice';
ax.fig = 'none';
tasks = {'B'};

% Choice + confidence: Fixed, Bayes, Lin, Quad
ah = show_data_or_fits('models', modelmaster([7 2 5 6]), 'root_datadir', datadir, 'depvars', {'resp','tf','Chat', 'g'}, 'axis', ax, ...
    'slices', {'c_s'}, 'tasks', tasks, 'symmetrify', false, 'errorbarwidth', 2.3, 'linewidth', 1.5,...
    'meanlinewidth', 1.5, 'gutter', [.01 .08], 'margins', [.1 .1 .1 .1], ...
    's_labels', [-9 -5 -2 0 2 5 9], 'show_legend', true, 'nPlotSamples', 2, 'nFakeGroupDatasets', 2,...
    'MCM', 'loopsis', 'MCM_size', .5, 'nRespSquares', 6, 'ref_model', 1, 'label_s_bin_centers', true, ...
    'sort_subjects', false, 'keep_subjects_sorted', true, 'show_subject_names', true, 'sort_model', 2, ...
    'reorder_subjects', reorder_subjects);

%     'sort_subjects', true, 'keep_subjects_sorted', true, 'show_subject_names', true, 'sort_model', 2);

h = gcf;
title(h.Children(5),'Quadratic')
xlabel(h.Children(5),'')
title(h.Children(6),'Linear')
xlabel(h.Children(6),'')
title(h.Children(7),'Bayesian')
xlabel(h.Children(7),'')
ylabel(h.Children(9),'Mean response')
xlabel(h.Children(9),'Orientation (°)')

hFitsAll = figure('Position', h.Position);
copyobj(h.Children, hFitsAll)

%% mean button press + fits and model comparison: Choice only
ax.col = 'model';
ax.row = 'slice';
ax.fig = 'none';
tasks = {'B'};

% Choice + confidence: Fixed, Bayes, Lin, Quad
ah = show_data_or_fits('models', modelmaster([13 8 11 12]), 'root_datadir', datadir, 'depvars', {'Chat'}, 'axis', ax, ...
    'slices', {'c_s'}, 'tasks', tasks, 'symmetrify', false, 'errorbarwidth', 2.3, 'linewidth', 1.5,...
    'meanlinewidth', 1.5, 'gutter', [.01 .08], 'margins', [.1 .1 .1 .1], ...
    's_labels', [-9 -5 -2 0 2 5 9], 'show_legend', true, 'nPlotSamples', 20, 'nFakeGroupDatasets', 1000,...
    'MCM', 'loopsis', 'MCM_size', .5, 'nRespSquares', 6, 'ref_model', 1, 'label_s_bin_centers', true, ...
    'sort_subjects', false, 'keep_subjects_sorted', true, 'show_subject_names', true, 'sort_model', 2, ...
    'reorder_subjects', reorder_subjects);

h = gcf;
title(h.Children(5),'Quadratic')
xlabel(h.Children(5),'')
title(h.Children(6),'Linear')
xlabel(h.Children(6),'')
title(h.Children(7),'Bayesian')
xlabel(h.Children(7),'')
ylabel(h.Children(9),'Prop. response "cat. 1"')
xlabel(h.Children(9),'Orientation (°)')

hFitsChoice = figure('Position', h.Position);
copyobj(h.Children, hFitsChoice)

%% SUPPLEMENT FIGURES
%% FIGURE S2
figS2H = figure('Name','FIG S2','Units','inches','Position',[1 10 11 8.5]);

% panel a: choice boundary illustration
axS2H(1:2) = copyobj(hCBT.Children, figS2H);

% panels b and c: model fits and model comp
srcAx = [1:7 9];
destAx = 3:11;
for i = 1:numel(srcAx)
    axS2H(destAx(i)) = copyobj(hFitsAll.Children(srcAx(i)), figS2H);
end

% set axis properties
for i = 1:numel(axS2H)
    axS2H(i).Units = 'inches';
    axS2H(i).FontSize = 14;
    axS2H(i).LineWidth = 1;
end

axS2H(1).Position = [1.35 7.5 1 1];
axS2H(2).Position = [1 6 3 2.5];

axS2H(3).Position = [8.5 0.2 2.25 2];
axS2H(4).Position = [6 0.2 2.25 2];
axS2H(5).Position = [3.5 0.2 2.25 2];
axS2H(6).Position = [1 0.2 2.25 2];

axS2H(7).Position = [8.5 3 2.25 2];
axS2H(8).Position = [6 3 2.25 2];
axS2H(9).Position = [3.5 3 2.25 2];
axS2H(10).Position = [1 3 2.25 2];

%% FIGURE S3
% Model recovery

%% FIGURE S4
figS4H = figure('Name','FIG S4','Units','inches','Position',[1 10 11 9.5]);

% panels a and b: model fits and model comp
srcAx = [1:7 9];
destAx = 1:8;
for i = 1:numel(srcAx)
    axS4H(destAx(i)) = copyobj(hFitsChoice.Children(srcAx(i)), figS4H);
end

% panel c: category decision boundaries for different models, one observer
axS4H(9:10) = copyobj(hCBExample.Children, figS4H);

% set axis properties
for i = 1:numel(axS4H)
    axS4H(i).Units = 'inches';
    axS4H(i).FontSize = 14;
    axS4H(i).LineWidth = 1;
end

axS4H(1).Position = [8.5 4.2 2.25 2];
axS4H(2).Position = [6 4.2 2.25 2];
axS4H(3).Position = [3.5 4.2 2.25 2];
axS4H(4).Position = [1 4.2 2.25 2];

axS4H(5).Position = [8.5 7 2.25 2];
axS4H(6).Position = [6 7 2.25 2];
axS4H(7).Position = [3.5 7 2.25 2];
axS4H(8).Position = [1 7 2.25 2];

axS4H(9).Position = [1.35 2.5 1 1];
axS4H(10).Position = [1 1 3 2.5];


%% grab params
model = 2; % Bayesian
paramName = 'logsigma_d';

paramIdx = find(strcmp(modelmaster(model).parameter_names, paramName));
nSubjects = numel(modelmaster(model).extracted);

vals = [];
for iS = 1:nSubjects
    vals(iS) = modelmaster(model).extracted(iS).best_params(paramIdx);
end

vals = exp(vals);

valsMedian = median(vals);
valsQ = prctile(vals,[25 75]);

