% meg_plotDecode.m

%% Setup
% load data_hdr.mat
classNames = {'V','H'};
p = [];

% fake data for testing
nTrials = 100;
classLabels = repmat([1 0],1,nTrials/2);
data = rand(200,157,nTrials);
data(:,:,classLabels==1) = data(:,:,classLabels==1)*2;

%% Decode
A = meg_decode(data, classLabels, classNames, p);

%% Unpack analysis structure
times = A.classTimes;
classAcc = A.classAcc;
classWeights = A.classWeights;

%% Plot
%% plot
xlims = targetWindow;
ylims = [30 100];

fH(1) = figure;
hold on
plot(times, classAcc, 'k')
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
            vals = classWeights(:,iTime,iC)';
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
            vals = squeeze(mean(abs(classWeights(:,tidx,iC,:)),2))';
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

vals = squeeze(mean(mean(abs(classWeights(:,tidx,:,:)),2),3))';
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
