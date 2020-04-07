function A = meg_decode(data, classLabels, classNames, p)

% function A = meg_decode(data, classLabels, classNames, p)
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
%   targetWindow: [start end] time window in which to decode
%   channels: channels to use for decoding, treated as indices
%   data_hdr: ft header for plotting topographies
%
% OUTPUTS
% A
%   analysis structure
% 
% Rachel Denison
% March 2020

%% Inputs
if nargin < 2 || isempty(data) || isempty(classLabels)
    error('data and classLabels are required inputs')
end
if nargin < 3 || isempty(classNames)
    classNames = {'class1','class2'};
end
if nargin < 4 || isempty(p)
    p = struct([]);
end

%% Unpack params structure
if isfield(p,'t')
    t = p.t;
else
    t = 1:size(data,1);
end
if ~isfield(p,'targetWindow') || isempty(p.targetWindow)
    targetWindow = t([1 end]); 
else
    targetWindow = p.targetWindow;
end
if ~isfield(p,'channels') || isempty(p.channels)
    channels = 1:size(data,2);
else
    channels = p.channels;
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
if size(classLabels,1)==1
    classLabels = classLabels'; % should be column vector
end

%% Organize data
dataInput = [];
for iC = 1:nC
    class = classes(iC);
    
    wC = classLabels==class; % select data in class
%     wNaN = isnan(squeeze(data(1,1,:))); % remove nan data
%     w = wC & ~wNaN;
    w = wC;
    
    dataInput{iC} = data(:,channels,w);
end

%% Decoding setup
getWeights = 1;

nt = 5; % 5 % average this many trials together to improve SNR
sp = 5; % 5 % sampling period
nReps = 10; % 10;
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);
svmopsNoCV = '-s 0 -t 0 -c 1 -q';
decodeAnalStr = sprintf('sp%d_nt%d_r%d', sp, nt, nReps);

%% Decoding
times = targetWindow(1):sp:targetWindow(2)-sp;

classAccNT = [];
for iRep = 1:nReps
    classAcc = [];
    vals1 = dataInput{1};
    vals2 = dataInput{2};

    % average trials
    if nt > 1
        vals1a = []; vals2a = [];
        n = size(vals1,3);
        trialsIdx = randperm(n);
        startTrials = 1:nt:n-nt;
        for iST = 1:numel(startTrials)
            trIdx = trialsIdx(startTrials(iST):startTrials(iST)+nt-1);
            vals1a(:,:,iST) = nanmean(vals1(:,:,trIdx),3);
            vals2a(:,:,iST) = nanmean(vals2(:,:,trIdx),3);
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
        
        % check
        if size(Xs,1) < kfold
            error('Size of data < kfold at time = %d', time)
        end
        
        % fit and cross validate classifier
        acc(iTime) = svmtrain(Y, Xs, svmops);
        
        % get the svm model, no cv
        if getWeights
            model(iTime) = svmtrain(Y, Xs, svmopsNoCV);
        else
            model = [];
        end
    end

    classAccNT(:,iRep) = acc;
    classModelNT{iRep} = model;
end

% trial average
classAcc = mean(classAccNT,2);

%% extract channel weights 
if getWeights
    classWeightsNT = [];
    for iTime = 1:numel(times)
        for iRep = 1:nReps
            model = classModelNT{iRep}(iTime);
            
            w = model.SVs' * model.sv_coef;
            b = -model.rho;
            if (model.Label(1) == -1)
                w = -w; b = -b;
            end
            classWeightsNT(iTime,:,iRep) = w;
        end
    end
    classWeights = mean(classWeightsNT,3);
else
    classWeights = [];
end

%% store results
A.classNames = classNames;
A.targetWindow = targetWindow;
A.decodingOps.channels = channels;
A.decodingOps.nTrialsAveraged = nt;
A.decodingOps.binSize = sp;
A.decodingOps.kfold = kfold;
A.decodingOps.nReps = nReps;
A.decodingOps.svmops = svmops;
A.decodingOps.analStr = decodeAnalStr;
A.classTimes = times;
% A.classAccNT = classAccNT;
% A.classModelNT = classModelNT;
% A.classWeightsNT = classWeightsNT;
A.classAcc = classAcc;
A.classWeights = classWeights;

