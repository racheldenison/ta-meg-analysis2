function [dataKT] = meg_subjectNames_RD2KT(dataRD)
% function [dataKT] = subjectNames_RD2KT(dataRD,subjectIdx)
% Inputs:
%   - dataRD: data matrix organized in RD subject order 

% Rearrange data matrix from RD subject odering to KT subject ordering 
% dataRD must be 1 x 10? 

% Subject oder RD 
% 1. R0817
% 2. R1187
% 3. R0983 
% 4. R0898
% 5. R1021
% 6. R1103
% 7. R0959
% 8. R1373
% 9. R1452
% 10. R1507

% Subject order KT 
% 1. R0817 
% 2. R0898 
% 3. R0959 
% 4. R0983
% 5. R1021
% 6. R1103 
% 7. R1187
% 8. R1373 
% 9. R1452
% 10. R1507

%% Check inputs
dims = size(dataRD); 
dimsProd = prod(dims); 

% dataRD_1D = reshape(dataRD,[dimsProd 1]); 

% dataKT = NaN(size(dataRD)); 

% testing reshaping 
% dataKT = NaN(dims);
% for i = 1:10
%     for iC = 1:2
%         for iT = 1:2
%             dataKT(iC,iT,i) = i*100 + iT*10 + iC;
%         end
%     end
% end
% dataKT_1 = reshape(dataKT,[dimsProd 1]);
% dataKT_2 = reshape(dataKT_1,[dims]);

%% 
subjectNamesKT = {'R0817',... 1
    'R0898',... 2
    'R0959',... 3
    'R0983',... 4
    'R1021',... 5
    'R1103',... 6
    'R1187',... 7
    'R1373',... 8 
    'R1452',... 9 
    'R1507',... 10 
    };

subjectNamesRD = {'R0817',... 1
    'R1187',... 2
    'R0983',... 3
    'R0898',... 4
    'R1021',... 5
    'R1103',... 6
    'R0959',... 7
    'R1373',... 8
    'R1452',... 9
    'R1507',... 10
    };

%% 
dataKT = []; 
for iS = 1:numel(subjectNamesRD)
    [~,Locb] = ismember(subjectNamesRD,subjectNamesKT{iS});
    sIdx = find(Locb); 
    dataKT(sIdx) = dataRD(iS); 
end



