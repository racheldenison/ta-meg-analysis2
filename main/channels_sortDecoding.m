function groupC = channels_sortDecoding
% Export channels by decoding weight 
% January 20, 2023

%% from Jiating (slack) 
% average top 50 weighted channels
% chsSorted = [ 41, 7, 111, 78, 36, 104, 92, 63, 54, 50, 20, 1, 16, 52, 8, 25, 60, 10, 26, 103, 15, 19, 34, 53, 32, 75, 44, 24, 2, 45, 37, 13, 82, 98, 65, 62, 29, 105, 110, 51, 81, 22, 73, 14, 59, 33, 55, 42, 28, 39];

% read text file of top 50 weighted channels x 20 sessions at 250 ms post
% T1
chFile = '/Users/kantian/Dropbox/Data/TA2/MEG/Group/channels/ch50_20ss.csv'; 
chMat = readtable(chFile);

%% 
clear groupC % group channels variable 
for i = 1:20 % sessions 
    groupC(i).channelsRanked = chsSorted; 
end

