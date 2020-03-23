function [A, fH] = meg_decode(data, classLabels, classNames, p)

% function [A, fH] = meg_decode(data, classLabels, classNames, p)
%
% INPUTS
% data
%   data matrix, time x channels x trials
%
% classLabels
%   vector of numerical class labels, trials x 1
%   NaN values are ignored.
% 
% classNames
%   cell array of class names (optional)
%
% p
%   params structure (optional)
%   t: time vector corresponding to data
%   targetWindow: time window in which to decode
%   channels: channels to use for decoding, treated as indices
%   data_hdr: ft header for plotting topographies
%
% OUTPUTS
% A
%   analysis structure
%
% fH
%   figure handles
% 
% Rachel Denison
% March 2020

%% Deal with inputs
if nargin < 4
    p = struct([]);
end
if nargin < 3
    classNames = {'class1','class2'};
end
if nargin < 2
    error('data and classLabels are required inputs')
end

%% Unpack params structure
if isfield(p,'t')
    t = p.t;
else
    t = 1:size(data,1);
end
if isfield(p,'targetWindow')
    targetWindow = p.targetWindow; 
else
    targetWindow = t([1 end]);
end
if isfield(p,'channels')
    channels = p.channels;
else
    channels = 1:size(data,2);
end
if isfield(p,'hdr')
    data_hdr = p.data_hdr; % data header for plotting topologies
else
    load data/data_hdr.mat
end

%% Setup
classes = unique(classLabels(~isnan(classLabels)));
nC = numel(classNames);

%% Checks
if numel(t)~=size(data,1)
    error('Length of t does not match first dimension of data')
end
if nC~=2
    error('This function currently handles only 2 classes')
end
if numel(classes)~=nC
    error('Number of classes in classLabels does not match classNames')
end

%% Organize data
dataInput = [];
for iC = 1:nC
    class = classes(iC);
    
    wC = classLabels==class; % select data in class
    wNaN = isnan(squeeze(data(1,1,:))); % remove nan data
    w = wC & ~wNaN;
    
    dataInput{iC} = data(:,channels,w);
end

%% Decoding setup
getWeights = 1;
syntheticTrials = 0;

nSynTrials = 100; % if constructing synthetic trials
nt = 5; % 5 % average this many trials together to improve SNR
sp = 5; % 5 % sampling period
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);
svmopsNoCV = '-s 0 -t 0 -c 1 -q';
decodeAnalStr = sprintf('sp%d_nt%d', sp, nt);

if syntheticTrials
    nReps = 1;
else
    nReps = nt;
end

%% Decoding
times = targetWindow(1):sp:targetWindow(2);

vals1 = dataInput{1};
vals2 = dataInput{2};

classAccNT = [];
for iRep = 1:nReps
    classAcc = [];

    % average trials
    if nt > 1
        vals1a = []; vals2a = [];
        n = size(vals1,3);
        if syntheticTrials
            nIdx = nSynTrials*nt;
            trialsIdx = [];
            for i = 1:ceil(nIdx/n)
                trialsIdx = [trialsIdx randperm(n)];
            end
            startTrials = 1:nt:nIdx;
        else
            trialsIdx = randperm(n);
            startTrials = 1:nt:n-nt; % n -> n-nt
        end
        for iST = 1:numel(startTrials)
            trIdx = trialsIdx(startTrials(iST):startTrials(iST)+nt-1);
            vals1a(:,:,iST) = mean(vals1(:,:,trIdx),3);
            vals2a(:,:,iST) = mean(vals2(:,:,trIdx),3);
        end
        vals1 = vals1a; vals2 = vals2a;
    end
    
    vals0 = cat(3, vals1, vals2);
    labels0 = [ones(size(vals1,3),1); zeros(size(vals2,3),1)];
    
    %% stratify
    nSamples = numel(labels0);
    foldSize = ceil(nSamples/kfold/2); % 2 classes
    stratIdx = [];
    for iFold = 1:kfold
        idx1 = (1:foldSize) + (iFold-1)*foldSize;
        idx2 = idx1 + nSamples/2;
        stratIdx = [stratIdx idx1 idx2];
    end
    stratIdxS = sort(stratIdx);
    r = stratIdxS(diff(stratIdxS)==0);
    ridx = [];
    for iR = 1:numel(r)
        ridx(iR) = find(stratIdx==r(iR),1,'last');
    end
    stratIdx(ridx) = [];
    if numel(stratIdx)>numel(labels0)
        stratIdx(numel(labels0)+1:end) = [];
    end
    
    vals = vals0(:,channels,stratIdx);
    labels = labels0(stratIdx);
    
    %% classify
    tic
    acc = [];
    for iTime = 1:numel(times)
        fprintf(' ')
        time = times(iTime);
        
        % classification data
        X = squeeze(mean(vals(find(t==time):find(t==time+sp-1),:,:),1))'; % average across time window
        Y = labels;
        
        % remove nan
        idx = isnan(X(:,1));
        X(idx,:) = [];
        Y(idx) = [];
        
        % scale data
        Xs = zscore(X);
        
        % fit and cross validate classifier
        acc(iTime) = svmtrain(Y, Xs, svmops);
        
        % get the svm model, no cv
        if getWeights
            model(iTime) = svmtrain(Y, Xs, svmopsNoCV);
        else
            model = [];
        end
    end
    toc
    
    classAcc(:,iC) = acc;
    classModel{iC} = model;
    
    % trial average
    classAccNT(:,:,iRep) = classAcc;
    classModelNT(:,iRep) = classModel';
end

%% extract channel weights 
if getWeights
    classWeights = [];
    for iC = 1:nC/2
        for iTime = 1:numel(times)
            for iRep = 1:nReps
                model = classModelNT{iC,iRep}(iTime);
                
                w = model.SVs' * model.sv_coef;
                b = -model.rho;
                if (model.Label(1) == -1)
                    w = -w; b = -b;
                end
                classWeights(:,iTime,iC,iRep) = w;
            end
        end
    end
else
    classWeights = [];
end

%% plot
xlims = targetWindow;
ylims = [30 100];

fH(1) = figure;
hold on
plot(times, mean(classAccNT,3),'LineWidth',1)
plot(times, mean(mean(classAccNT,3),2), 'k')
plot(xlims,[50 50],'k')
xlim(xlims)
ylim(ylims)
xlabel('Time (ms)')
ylabel('Classification accuracy (%)')
legend(classNames)

% if saveFigs
%     rd_saveAllFigs(gcf, {sprintf('%s_%s',figName{1},decodeAnalStr)}, 'plot', figDir)
% end

%% topo weights movie T1 and T2
if getWeights
    clims = [-2 2];
    
    fH(2) = figure('Position',[250 850 950 450]);
    for iTime = 1:numel(times)
        for iC = 1:nC/2
            vals = squeeze(mean(classWeights(:,iTime,iC,:),4))';
            subplot(1,nC/2,iC)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(classNames{iC})
        end
        rd_supertitle2(sprintf('t = %d', times(iTime)))
        pause(0.2)
    end
end

%% topo weights for specific time intervals
twins = {[110 140], [140 230], [230 280], [280 315], [110 315]};

if getWeights
    clims = [0 1.5];

    for iTW = 1:numel(twins)
        twin = twins{iTW};
        tidx = find(times==twin(1)):find(times==twin(2));
        
        figure('Position',[250 850 950 450]);
        for iC = 1:nC/2
            vals = squeeze(mean(mean(abs(classWeights(:,tidx,iC,:)),4),2))';
            subplot(1,nC/2,iC)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(classNames{iC})
        end
        rd_supertitle2(sprintf('%d-%d ms', twin(1), twin(2)))
        
%         if saveFigs
%             rd_saveAllFigs(gcf, ...
%                 {sprintf('%s_%s_%d-%dms','svmWeights',decodeAnalStr, twin(1), twin(2))}, 'map', figDir)
%         end
    end
end

%% mean across longest interval, reps, and orientations
twin = [110 315];
tidx = find(times==twin(1)):find(times==twin(2));
nTopChannels = 10;

vals = squeeze(mean(mean(mean(abs(classWeights(:,tidx,:,:)),4),2),3))';
[sortedVals, idx] = sort(vals,'descend');
topChannels = idx(1:nTopChannels);

fH(3) = figure('Position',[250 850 950 450]);
subplot(1,3,1)
histogram(vals)
xlabel('Unsigned SVM weight')
ylabel('Count')
subplot(1,3,2)
ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
title('Unsigned SVM weights')
subplot(1,3,3)
ssm_plotOnMesh(double(vals>=sortedVals(nTopChannels)), '', [], data_hdr, '2d');
title(['Channels ' sprintf('%d ',topChannels)])
rd_supertitle2(sprintf('%d-%d ms', twin(1), twin(2)))

% if saveFigs
%     rd_saveAllFigs(gcf, ...
%         {sprintf('%s_%s_%d-%dms_top%dCh','svmWeights',decodeAnalStr, twin(1), twin(2), nTopChannels)},...
%         'map', figDir)
% end

%% store results
A.classNames = classNames;
A.targetWindows = targetWindow;
A.decodingOps.channels = channels;
A.decodingOps.nTrialsAveraged = nt;
A.decodingOps.binSize = sp;
A.decodingOps.kfold = kfold;
A.decodingOps.svmops = svmops;
A.classTimes = times;
A.classAcc = classAcc;
A.classModel = classModel;
A.classWeights = classWeights;

% %% save analysis
% if saveAnalysis
%     save(sprintf('%s_%s_%s.mat',analysisFileName,analStr,decodeAnalStr), 'A')
% end
