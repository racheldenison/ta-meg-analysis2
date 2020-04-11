function [D, I] = meg_slicer(data, cond, condNames, levelNames)

% function [D, I] = meg_slicer(data, cond, condNames, levelNames)
%
% INPUTS
% data
%   data matrix, time x channels x trials
%
% cond
%   condition matrix, trials x conditions
%   Each column is the condition associated with one slice. The value
%   gives the (categorical) level of that condition. NaN values are
%   excluded.
%
% condNames
%   condition names cell array, 1 x conditions
%
% levelNames
%   level names cell array, 1 x conditions
%   Each entry is another cell array giving the level names for that
%   condition.
%
% See "args" and "set up data" cells for examples.
%
% OUTPUTS
% D
%   slice data structure
%   Each field is a condition combination containing a data matrix (time x 
%   channels x trials) for the trials in that group.
%
% I
%   trial indices for each condition combination
%
%
% Rachel Denison
% January 2020

%% args
if nargin < 1
    generateData = 1;
else
    generateData = 0;
end

if nargin < 2 && ~generateData
    error('Must input cond with data')
end

if nargin < 3
    condNames = {'cue','T1','T2'};
%     condNames = {'cue'};
end

if nargin < 4
    levelNames = {{'cueT1','cueT2','cueN'},{'t1v','t1h'},{'t2v','t2h'}};
%     levelNames = {{'cueT1','cueT2'}};
end

%% check inputs
% check data
if ~generateData
    if isempty(data)
        D = [];
        I = [];
        return;
    end
    sz = size(data);
    if numel(sz)~=3
        error('data is expected to be 3-dimensional, with trials as the last dimension')
    end
end

% check levelNames
if numel(levelNames)~=numel(condNames)
    error('Must be equal number of condNames and levelNames')
end
if ~iscell(levelNames{1})
    error('levelNames should be a cell array of cell arrays')
end

%% set up data
if generateData
    % generate fake data
    nT = 1000;
    nCh = 157;
    nTrialsPerCond = 100;
    
    nLevels = [];
    for iSlice = 1:numel(condNames)
        nLevels(iSlice) = numel(levelNames{iSlice});
    end

    cond0 = fullfact([nLevels nTrialsPerCond]);
    cond = cond0(:,1:end-1);
    nTrials = size(cond,1);
    
%     data = rand(nT,nCh,nTrials);
    % fake data where the value is the trial number
    data0 = shiftdim(1:nTrials,-1);
    data = repmat(data0,nT,nCh);
end

%% set up slices
nSlices = numel(condNames);

w = []; % logical matrix for each slice, nTrials x nLevels
nLevels = [];
for iSlice = 1:nSlices
    sliceCond = cond(:,iSlice);
    levels = unique(sliceCond(~isnan(sliceCond)));
    nLevels(iSlice) = numel(levels);
    
    if nLevels ~= numel(levelNames{iSlice})
        error('The number of unique values in each column of cond should match the number of level names.')
    end
    
    for iLevel = 1:nLevels(iSlice)
        w{iSlice}(:,iLevel) = sliceCond==levels(iLevel);
    end
end

%% slice data
combos = fullfact(nLevels); % each column is a slice, and the value is the level of that slice
nCombos = size(combos, 1);

for iCombo = 1:nCombos
    combo = combos(iCombo,:);
    
    wSlices = [];
    comboName = '';
    for iSlice = 1:nSlices
        level = combo(iSlice);
        wSlices(:,iSlice) = w{iSlice}(:,level);
        
        levelName = levelNames{iSlice}{level};
        comboName = strcat(comboName,levelName,'_');
    end
    comboName(end) = [];
    
    wCombo = sum(wSlices,2)==nSlices; % all conditions are met
        
    D.(comboName) = data(:,:,wCombo);
    I.(comboName) = find(wCombo);
end
     
