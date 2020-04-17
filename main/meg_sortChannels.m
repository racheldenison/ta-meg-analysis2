function [C, fH, figNames] = meg_sortChannels(expt, sessionDir, user, sortType)

%% setup
%% inputs
if nargin == 0
    % for testing
    expt = 'TA2';
    sessionDir = 'R1187_20181119';
    user = 'mcq'; 
    sortType = 'classweights'; % alpha, threshold
end

% directories
exptDir = meg_pathToTAMEG(expt, user);

% params
paramType = 'Analysis';
p = meg_params(sprintf('%s_%s', expt, paramType));

% topo header
load vis/data_hdr.mat

%% sort channels
%%% RD: Add Karen's other sortTypes here from meg_selectChannels
% make rank-ordered list of channels, channelsRanked
switch sortType
    case 'peakprom'
        
    case 'classweights'
        %%% RD: Make separate function
        twin = [0 400];
        
        matDir = sprintf('%s/%s/mat', exptDir, sessionDir);
        dataFile = dir(sprintf('%s/classAcc*varyNChannels*.mat', matDir));
        load(sprintf('%s/%s', matDir, dataFile.name))

        A = allA{1};
        nT = numel(A);
        c = 'all';
        t = A(1).(c).classTimes - p.eventTimes(2);
        tidx = find(t==twin(1)):find(t==twin(2));
        
        classWeights = [];
        for iT = 1:nT
            classWeights(:,:,iT) = A(iT).(c).classWeights(tidx,:,:);
        end
        
        absClassWeight = abs(mean(mean(classWeights,3),1))';
        [weightsSorted, sortIdx] = sort(absClassWeight,1,'descend');
        
        channelsRanked = sortIdx;
        
        C.sortType = sortType;
        C.sessionDir = sessionDir;
        C.fileName = dataFile.name;
        C.twin = twin;
        C.classWeights = classWeights;
        C.absClassWeight = absClassWeight;
        C.weightsSorted = weightsSorted;
        C.channelsRanked = channelsRanked;
        
    case 'alpha'
        
    otherwise
        error('sortType not recognized')
end

%% plot
switch sortType
    case 'classweights'
        chVal = absClassWeight';
end
chRankVal = linspace(1,0,numel(channelsRanked));
[Y, idx] = sort(channelsRanked);

fH = figure;
subplot(1,2,1)
ssm_plotOnMesh(chVal, '', [], data_hdr, '2d');
colorbar
title('abs class weight')
subplot(1,2,2)
ssm_plotOnMesh(chRankVal(idx), '', [], data_hdr, '2d');
colorbar
title('by rank')

figNames = {sprintf('map_channelsRanked_%s', sortType)};

