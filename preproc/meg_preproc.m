function preprocFileName = meg_preproc(filename, figDir, badChannels)

%% Setup
% desk Rachel
% filename = '/Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/R0817_TAPilot_8.20.14.sqd';
% filename = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/preproc/R0890_TAPilot_8.06.14_run01.sqd';
% figDir = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/Runs/figures';

% desk Karen
% filename = '/Local/Users/kantian/Dropbox/Data/TAPilot/MEG/R0817_20140820/R0817_TAPilot_8.20.14.sqd';
% filename = '/Users/kantian/Dropbox/Data/TA2/MEG/R1507_20190627/R1507_TA2_6.27.19_run01.sqd';
% figDir = 'Users/kantian/Dropbox/Data/TA2/MEG/R1507_20190627/preproc/figures';

% remember, these channel numbers use zero indexing
megChannels = 0:156;
refChannels = 157:159;
triggerChannels = 160:167;
eyeChannels = 176:177;
photodiodeChannel = 191;

% badChannels = [];
if nargin < 3
    badChannels = [];
end

%%  preproc options

exptShortName = 'TANoise_Preproc'; % TA2_Preproc, TANoise_Preproc, Cupcake
p = meg_params(exptShortName); 

environmentalDenoise = 1; % e 
removeBadChannels = 1; % b dead and outlier sd channels 
interpolate = 1; % i interpolate bad channels 

TSPCA = 0; % t time shift pca 

hpfilter = 0; % f high pass filter 

components = 1; % c, pca/ica
rejectPC = 1; % auto reject 1st pc?
rejectIC = 0; % auto reject ic?

%%
% high pass filter options
Fsample = 1000;
Fhp = 0.1; % 1, 0.1 high pass frequency
N = []; %16500; % 8250 % filter order, auto calculate if unspecified 
type = 'firws';
direc = 'onepass-zerophase';

% trial definition (for pca/ica)
trialDef.trialFunHandle = @mytrialfun_all;
trialDef.prestim = p.prestim; % preproc timing window (different from analysis window) 
trialDef.poststim = p.poststim; 
trialDef.trig = p.trialDefTrig; 
trialDef.nTrigsExpected = [];
trialDef.threshold = 2.5; 

plotFigs = 1;
saveFigs = 1;

analStr = [];

deleteStepFiles = 1; 

%% Get the MEG data
% data is time x channels
[data, info] = sqdread(filename);

dataOriginal = data; % unperturbed copy for plots 

Fs = info.SampleRate;
t = 0:1/Fs:size(data,1)/Fs-1/Fs;

%% Set aside data from special channels
% add 1 to adjust for zero-indexing

refData = data(:,refChannels+1);
trigData = data(:,triggerChannels+1);
eyeData = data(:,eyeChannels+1);
pdData = data(:,photodiodeChannel+1);

%% Look at special channels
if plotFigs
    figure
    subplot(4,1,1)
    plot(t,refData)
    title(['reference channels' num2str(refChannels)])
    subplot(4,1,2)
    plot(t,trigData)
    legend(num2str(triggerChannels'))
    title(['trigger channels' num2str(triggerChannels)])
    subplot(4,1,3)
    plot(t,eyeData)
    title(['eye channels' num2str(eyeChannels)])
    xlabel('time (s)')
    subplot(4,1,4)
    plot(t,pdData)
    title(['photodiode channel' num2str(photodiodeChannel)])
    xlabel('time (s)')
end

%% Denoise using reference channels (new meg-utils)
if environmentalDenoise 
    % See also LSdenoise.m
    analStr = [analStr 'e'];
    
    opts = [];
    opts.ref = refChannels+1;
    opts.meg = megChannels+1;
    
    % convert to time x trials x channels
    data = permute(data,[1 3 2]);
    data = meg_environmentalDenoising(data, opts);
    data = permute(data,[1 3 2]); % convert back
    
    if plotFigs
        figure
        hold on
        subplot 211
        plot(dataOriginal(:,megChannels+1))
        title('unprocessed')
        subplot 212
        plot(data(:,megChannels+1)) 
        title('environmental denoise')
    end
    
%     % delete dead / interp file if exist
%     if deleteStepFiles
%         if exist(badFile, 'file')
%             delete(badFile);
%         elseif exist(interpFile, 'file')
%             delete(interpFile)
%         end
%     end
%     

    % write and save sqd    
    edFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    
    if exist(edFile,'file')
        error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
    else
        sqdwrite(filename, edFile, 'data', data);
    end
    
    dataset = edFile;
end

%% Find bad channels
if removeBadChannels
    analStr = [analStr 'b'];
    
    % high or low variance across the entire time series
    outlierSDChannels = meg_find_bad_channels(permute(data(:,megChannels+1),[1 3 2]));
    
    % dead or saturating channels for all or portions of the time series
    deadChannels = checkForDeadChannels(filename)+1;
    
    % aggregate the bad channels
    badChannels = unique([badChannels outlierSDChannels' deadChannels]);
    nBad = numel(badChannels);
    
    % plot the time series for the bad channels
    if plotFigs
        figure
        for iBad = 1:nBad
            subplot(nBad,1,iBad)
            chan = badChannels(iBad);
            plot(t, data(:,chan))
            title(sprintf('channel %d', chan))
        end
        xlabel('time (s)')
    end
    
    % zero out bad channels
    data(:,badChannels) = 0;

    % write and save sqd    
    badFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    
    if exist(badFile,'file')
        error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
    else
        sqdwrite(filename, badFile, 'data', data);
    end
    
    if deleteStepFiles
        delete(edFile);
    end
    
    dataset = badFile;
end

%% Interpolate to replace bad channels 
if interpolate 
    % create dummy ft_data from the original data and update it
%     dataset = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    analStr = [analStr 'i'];
    
    % read data into ft structure in continuous mode by initializing cfg with
    % only the dataset
    % ft data is channels x time 
    cfg = [];
    cfg.dataset = dataset;
    ft_data = ft_preprocessing(cfg); % this preserves the data, but scales it by like 10^-13
    ft_data.trial{1} = data'; % so just replace data with the original data
    
    cfg = [];
    cfg.badchannel = ft_channelselection(badChannels, ft_data.label);
    cfg.method = 'spline';
    % cfg.neighbours = neighbours; % only needed for nearest and average, not for spline
    ft_data = ft_channelrepair(cfg, ft_data);
    
    % compare before and after interpolation
    if plotFigs
        figure
        sampleChannels = badChannels;
        for iCh = 1:numel(sampleChannels)
            chan = sampleChannels(iCh);
            subplot(numel(sampleChannels),1,iCh)
            hold on
            plot(t, data(:,chan))
            plot(t, ft_data.trial{1}(chan,:), 'r')
            xlabel('time (s)')
            title(sprintf('channel %d', chan))
        end
        legend('before interpolation','after interpolation')
    end
   
    % save the interpolated data
    data(:,1:numel(megChannels)) = ft_data.trial{1}';
    
    % write interp file 
    interpFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    if exist(interpFile,'file')
        error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
    else
        sqdwrite(filename, interpFile, 'data', data);
    end
    
    % delete dead file if exist
    if exist(badFile, 'file')
        delete(badFile);
    end
    
    dataset = interpFile; 
    badChannels = []; 
end

%% Save data preprocessed up to this point
% preFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
% 
% if exist(preFile,'file')
%     error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
% else
%     sqdwrite(filename, preFile, 'data', data);
% end
% 
% dataset = preFile;

%% Time-shift PCA for environmental denoising
% http://www.isr.umd.edu/Labs/CSSL/simonlab/Denoising.html
% http://lumiere.ens.fr/Audition/adc/meg/
% "Noise fields measured by reference magnetometers are optimally filtered 
% and subtracted from brain channels. The filters (one per reference/brain 
% sensor pair) are obtained by delaying the reference signals, 
% orthogonalizing them to obtain a basis, projecting the brain sensors onto 
% the noise-derived basis, and removing the projections to obtain clean 
% data." (DOI: 10.1016/j.jneumeth.2007.06.003)
if TSPCA
    sizeOfBlocks = 20000;
    shifts = -100:100;
    sourceFile = preFile;
    
    analStr = [analStr 't'];
    tspcaFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    
    % run sqd denoise
    % this writes the tspca sqd file
    if exist(tspcaFile,'file')
        error('tspcaFile already exists ... will not overwrite. note that continuing the script will delete this file.')
    else
        fprintf('Running sqdDenoise\n');
        sqdDenoise(sizeOfBlocks, shifts, 0, sourceFile, [], 'no', ...
            p.channelForSaturatingChannels, 'yes', tspcaFile); % do not zero saturated 
        % sqdDenoise(sizeOfBlocks, shifts, 0, sourceFile, badChannels-1, 'no', ...
            % precueChannel, 'yes', tspcaFile); % do not zero saturated 
    end
    
    % plot comparison with basic environmental denoising
    if plotFigs
        datats = sqdread(tspcaFile);
        figure
        sampleChannels = [1 14 badChannels];
        for iCh = 1:numel(sampleChannels)
            chan = sampleChannels(iCh);
            subplot(numel(sampleChannels),1,iCh)
            hold on
            plot(t, data(:,chan))
            plot(t, datats(:,chan), 'r')
            xlabel('time (s)')
            title(sprintf('channel %d', chan))
        end
        legend('LS environmental denoise','+time-shift PCA environmental denoise')
    end
    
    % replace data with datats
    data = datats;
    % clear datats;
    
    % remove just-created sourceFile
    if deleteStepFiles
        delete(sourceFile);
    end
   
    dataset = tspcaFile;
end

%% High pass filter
if hpfilter
    analStr = [analStr 'f'];
    
    % data should be channels x time
    data = ft_preproc_highpassfilter(data', Fsample, Fhp, N, type, direc); 
    data = data';
    
    if plotFigs
        figure
        hold on
        subplot 211
        plot(dataOriginal(:,megChannels+1))
        title('unprocessed')
        subplot 212
        plot(data(:,megChannels+1))
        title('high pass filter')
    end

    filterFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    
    if exist(filterFile,'file')
        error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
    else
        sqdwrite(filename, filterFile, 'data', data);
    end
    
    if deleteStepFiles
        if TSPCA
            if exist(tspcaFile, 'file')
                delete(tspcaFile);
            end
        elseif interpolate 
            if exist(interpFile, 'file')
                delete(interpFile);
            end
        end
    end

    dataset = filterFile;
end

%% PCA/ICA
if components
    [ft_cleandata, WorkFlow, pcStr, ft_PCA] = meg_pca_ica(dataset, [], trialDef, rejectPC, rejectIC);
    analStr = [analStr pcStr];
    data(:,1:157) = ft_cleandata.trial{1}'./1e-15; 
    
    % save pca/ica data 
    pcaFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    if exist(pcaFile,'file')
        error('%s_%s.sqd already exists ... will not overwrite. exiting.', filename(1:end-4), analStr)
    else
        sqdwrite(filename, pcaFile, 'data', data);
    end
    
    clear ft_cleandata
    dataset = pcaFile; 
    
    if deleteStepFiles
        if exist(filterFile, 'file')
            delete(filterFile);
        end
    end
    
end

%% finally, check the triggers
meg_checkTriggers(dataset); % filename

%% return preproc file name
preprocFileName = sprintf('%s_%s.sqd', filename(1:end-4), analStr);

%% figure concat time series, channel avg trial time series, fft
img = meg_plotTime(preprocFileName, analStr, 0, 0, p); 

%% save figs
if saveFigs
    runStr = rd_getTag(filename);
    figSubDir = sprintf('%s/%s/%s', figDir, analStr, runStr);
    if ~exist(figSubDir,'dir')
        mkdir(figSubDir)
    end
    f = sort(findobj('Type','figure'));
    for iF = 1:numel(f)
        if isnumeric(f)
            figNames{iF} = num2str(f(iF));
        else
            figNames{iF} = num2str(f(iF).Number);
        end
    end
    rd_saveAllFigs(f,figNames,[],figSubDir);
    close all
end

