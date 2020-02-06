function data = meg_getData(sqdfile,p) 
% function D = MEG_GETDATA(filename)
%
% INPUT
%   filename: string path to sqd file 
%   p: structure of expt parameters 
% OUTPUT
%   data: time x channels x trials matrix 
% 
% Karen Tian 
% January 2020 

% dat = ft_read_data(sqdfile);
% hdr = ft_read_header(sqdfile);

tic
cfg                     = [];
cfg.dataset             = sqdfile;
cfg.trialdef.prestim    = p.prestim; %0; %0.5; % sec
cfg.trialdef.poststim   = p.poststim; %2.3; %3.1;
cfg.trialdef.trig       = [168,167]; %[168,167]; %[161:164,167]; %168 = precue, 167=blank
threshold = 2.5;
[trl,Events]            = mytrialfun_all(cfg,threshold,[]);

prep_data          = ft_preprocessing(struct('dataset',sqdfile,...
    'channel','MEG','continuous','yes','trl',trl));
toc

%% read and concatenate data 
tic
input = prep_data; 

vals = [];
for iTrial = 1:length(input.trial) % 580 trials
        vals2 = input.trial{1,iTrial}; % go through trial by trial cells 
        vals = [vals vals2]; % concat horizontally into trial x channel matrix 
end
dataConcat = vals; % use this as unfiltered concatenated data matrix chans x time 
 
nTrials = size(input.trial,2);
nChannels = size(input.trial{1},1); 
nTime = size(input.trial{1},2); 
toc
%% reshape 

% reshape concat time series into trial time series 
toReshape = dataConcat'; 
reshape1 = toReshape;
reshape2 = reshape(reshape1,nTime,nTrials,nChannels); 
reshape3 = permute(reshape2,[1 3 2]); 
reshape4 = mean(reshape3,3); 
data = reshape3; % time x channels x trials 

end