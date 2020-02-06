function meg_multiplot(vals, layout, cmap, highlightChannels)
%
% function meg_topoplot(vals, [layout], [cmap], [highlightChannels])
%
% INPUTS
% vals is an N x t matrix one value per channel, where N is number of channels
% layout is an ft layout (see ft_prepare_layout). will generate layout from
% custom file if not included.
% cmap is a colormap of any size. default is parula, which we load from
% saved mat file.
% highlightChannels is an Nx1 vector of channel numbers to be highlighted
%
% NOTE although layout and cmap are optional, this will run faster if you
% give them as arguments, because otherwise files must be loaded and layout
% must be prepared.

%% inputs
if nargin<2 || isempty(layout)
    load data_hdr.mat
    cfg = [];
    layout = ft_prepare_layout(cfg, data_hdr);
end
if nargin<3 || isempty(cmap)
    % load parula
    cmap = parula;
end
if nargin<4 || isempty(highlightChannels)
    highlight = 0;
else
    highlight = 1;
end

nChannels = 157;

%% cfg
cfg = [];
cfg.style = 'straight'; % no contours
cfg.interpolation = 'nearest'; % no smooth interpolation
if highlight
    cfg.highlight = 'label';
    cfg.highlightchannel = highlightChannels;
    cfg.highlightsize = 20; 
    cfg.highlightcolor = [1 0 0];  
    cfg.showlabels = 'yes'; 
    % cfg.marker = 'numbers';
else
    cfg.showlabels = 'yes'; 
    cfg.marker = 'numbers';
end
cfg.layout = layout;
cfg.colormap = cmap;
cfg.comment = 'no'; % no info about date, time, data range

%% data
p = meg_params('TA2');
data = [];
data.label = layout.label(1:nChannels);
data.avg = vals;
data.time = p.tstart:p.tstop;
data.dimord = 'chan_time';

%% plot
ft_multiplotER(cfg, data);

