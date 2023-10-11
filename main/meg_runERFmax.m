% Find ERF max 
% using top decoding channels 

%% Settings (to change) 
% nChOI = 5; % number of channels to select 

%% Settings (more fixed) 
% slice = 'cue'; 
% expt = 'TA2'; 
% user = 'karen'; 

%% load group data 
% data structure named D4 
% fields for each attention condition (cueT1, cueT2, cueN) 
% 4401 time x 157 channels x 20 sessions 
load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/mat/D4.mat')
fields = fieldnames(D4); 

%% Load channels by decoding weights 
% Align subject IDs 
chsJZ = load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/channels/ch50_20ss.csv'); % from Jiating, using different subject order

sessionIDsJZ = {'R0817_20181120', 'R0817_20190625',...
     'R1187_20181119', 'R1187_20190703',...
     'R0959_20181128', 'R0959_20190703',...
     'R1373_20181128', 'R1373_20190708',...
     'R0983_20190722', 'R0983_20190723',...
     'R1507_20190621', 'R1507_20190627',...
     'R0898_20190723', 'R0898_20190724',...
     'R1452_20181119', 'R1452_20190711',...
     'R1547_20190729', 'R1547_20190730',...
     'R1103_20181121', 'R1103_20190710'};

sessionIDsKT =  {'R0817_20181120', 'R0817_20190625',...
            'R0898_20190723', 'R0898_20190724',...
            'R0959_20181128', 'R0959_20190703',...
            'R0983_20190722', 'R0983_20190723',...
            'R1103_20181121', 'R1103_20190710',...
            'R1187_20181119', 'R1187_20190703',...
            'R1373_20181128', 'R1373_20190708',...
            'R1452_20181119', 'R1452_20190711',...
            'R1507_20190621', 'R1507_20190627',...
            'R1547_20190729', 'R1547_20190730'}; 

% Match index of session IDs 
for i = 1:numel(sessionIDsJZ)
    idKT = sessionIDsKT{i};
    index = find(contains(sessionIDsJZ,idKT));
    chs(i,:) = chsJZ(index,:)+1;  % correction for 0 index 
end

% Put channels in group channels structure groupC
groupCRef = load('/Users/kantian/Dropbox/Data/TA2/MEG/Group/mat/groupC.mat'); % load old groupC for channel Prominence 

clear groupC
for iS = 1:size(chs,1)
    groupC(iS).sessionDir = sessionIDsKT{i};
    groupC(iS).selectionMethod = 'decodingWeights';
    groupC(iS).channelsRanked = chs(iS,:);
    for iC = 1:size(chs,2)
        ch = chs(iS,iC); 
        groupC(iS).channelDirection(iC) = groupCRef.groupC(iS).chPromDir(find(groupCRef.groupC(iS).sortChByProm==ch));
    end
end

%% 
nChOI = 1:50; 
expt = 'TA2';
user = 'karen'; 

[D,A] = meg_peakAnalysis(D4,nChOI,groupC,expt,user); 

%% ANOVA quick 
count = 1; 
clear g 
targetNames = {'T1','T2'}; 
targets = 1:2; 
s = 1:10; 
att = 1:3; 
for iF = 1:numel(fields)
    for iT = 1:2
        for iS = 1:2:20
            y(count:count+1) = A.(targetNames{iT}).(fields{iF}).maxSession(iS:iS+1); 
            g.att(count:count+1) = [att(iF) att(iF)]; 
            g.subject(count:count+1) = [s((iS+1)/2) s((iS+1)/2)]; 
            g.target(count:count+1) = [targets(iT) targets(iT)]; 
            count = count+2; 
        end
    end
end

[p tbl stats] = anovan(y,{g.att,g.subject,g.target},'model','interaction','varnames',{'att','subject','target'});

yT1cT1 = y(g.target==1 & g.att==1); 
yT1cT2 = y(g.target==1 & g.att==2); 
[h.yT,p,ci,stats] = ttest2(yT1cT1,yT1cT2);

yT2cT1 = y(g.target==2 & g.att==1); 
yT2cT2 = y(g.target==2 & g.att==2); 
[h,p,ci,stats] = ttest2(yT2cT1,yT2cT2);






