function [A, fH, figNames] = meg_plotDecode(D, I, p, classLabels, classNames)

% function [A, fH, figNames] = meg_plotDecode(D, I, p, classLabels, classNames)
%
% ADD doc
%
% Rachel Denison
% March 2020

%% Inputs
if nargin < 1
    D = [];
end
if nargin < 2
    I = [];
end
if nargin < 2
    p = [];
end
if nargin < 3
    classLabels = [];
end
if nargin < 4 
    classNames = [];
end

%% Setup
if isfield(p,'hdr')
    data_hdr = hdr;
else
    load data_hdr.mat
end

%% Fake data if requested
if isempty(D)
    fprintf('\nGenerating fake data ...\n')
    nT = 200;
    nCh = 157;
    nTrials = 500;
    classLabels = repmat([1 0],1,nTrials/2);
    data = rand(nT,nCh,nTrials);
    data(:,:,classLabels==1) = data(:,:,classLabels==1)*1.1;
    data = data + randn([nT nCh nTrials])*1;
    D.fake = data;
    I.fake = 1:nTrials;
end

%% Get conditions
condNames = fields(D);
nCond = numel(condNames);

%% Decode
for iC = 1:nCond
    c = condNames{iC};
    data = D.(c);
    
    if isempty(I)
        idx = 1:size(data,3);
    else
        idx = I.(c);
    end

    A.(c) = meg_decode(data, classLabels(idx), classNames, p);
end

%% Unpack analysis structure
% same for all conds
classNames = A.(c).classNames;
times = A.(c).classTimes;
targetWindow = A.(c).targetWindow;
decodeAnalStr = A.(c).decodingOps.analStr;

classAcc = [];
classWeights = [];
for iC = 1:nCond
    c = condNames{iC};
    classAcc(:,iC) = A.(c).classAcc;
    classWeights(:,:,iC) = A.(c).classWeights;
end

%% Plot
fH = [];
figNames = {};
decodeTitle = sprintf('%s vs. %s', classNames{1},classNames{2});

%% plot
xlims = targetWindow;
ylims = [30 100];

fH(1) = figure;
hold on
plot(times, classAcc)
plot(xlims,[50 50],'k')
xlim(xlims)
ylim(ylims)
xlabel('Time (ms)')
ylabel('Classification accuracy (%)')
legend(condNames)
title(decodeTitle)

figNames{1} = sprintf('plot_%s_%s', 'classAcc', decodeAnalStr);

%% topo weights movie T1 and T2
plotMovie = 0;
if ~isempty(classWeights) && plotMovie
    clims = [-2 2];
    
    figure
    for iTime = 1:numel(times)
        subplot(1,1,1)
        vals = mean(classWeights(iTime,:,:),3);
        ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
        set(gca,'CLim',clims)
        colorbar
        rd_supertitle2(sprintf('%s, t = %d', decodeTitle, times(iTime)))
        pause(0.2)
    end
end

%% topo weights for specific time intervals
% twins = {[1000 1100], [1100 1200]};
twins = {targetWindow};

if ~isempty(classWeights)
    clims = [0 0.5];
    
    for iTW = 1:numel(twins)
        twin = twins{iTW};
        [~, tidx1] = min(abs(times-twin(1)));
        [~, tidx2] = min(abs(times-twin(2)));
        tidx = tidx1:tidx2;
        
        fH(end+1) = figure;
        vals = mean(mean(abs(classWeights(tidx,:,:)),1),3);
        ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
        set(gca,'CLim',clims)
        colorbar
        rd_supertitle2(sprintf('%s, %d-%d ms', decodeTitle, times(tidx1), times(tidx2)))
        
        figNames{end+1} = sprintf('map_svmWeights_%s_%d-%dms', decodeAnalStr, times(tidx1), times(tidx2));
    end
end

%% Save analysis
% if saveAnalysis
%     save(sprintf('%s_%s_%s.mat', analysisFileName, analStr, decodeAnalStr), 'A')
% end
