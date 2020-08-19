function [sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt)

% function [sessionNames,subjectNames,ITPCsubject,ITPCsession] = meg_sessions(expt)
%
% INPUT
% expt
%   experiment name, 'TANoise' or 'TA2'
%
% OUTPUT
% sessionNames
%   cell array of session names
% subjectNames
%   cell array of subject names 
% ITPCDir
%   array of subject ITPC evoked response direction (1 increase, -1 decrease)
% sessionITPCDir 

switch expt
    case 'TANoise'
        sessionNames = {'R0817_20171212','R0817_20171213',...
            'R0898_20180112','R0898_20180116',...
            'R0959_20180219','R0959_20180306',...
            'R0983_20180111','R0983_20180112',...
            'R1021_20180208','R1021_20180212',...
            'R1103_20180213','R1103_20180215',...
            'R1187_20180105','R1187_20180108',...
        	'R1373_20190723','R1373_20190725',...
            'R1452_20190717','R1452_20190718',...
            'R1507_20190702','R1507_20190705'}; % N=10 x 2 sessions TANoise
         subjectNames = {'R0817',...
            'R0898',...
            'R0959',...
            'R0983',...
            'R1021',...
            'R1103',...
            'R1187',...
        	'R1373',...
            'R1452',...
            'R1507'};
        ITPCsubject = [1,... % 'R0817'
            1,... % 'R0898'
            -1,...% 'R0959'
            -1,...% 'R0983'
            0,...% 'R1021'
            1,...% 'R1103'
            -1,...% 'R1187'
        	1,...% 'R1373'
            -1,...% 'R1452'
            -1% 'R1507'};
            ]; 
        
        subject = 1;
        for i = 1:2:numel(sessionNames)
            val = ITPCsubject(subject);
            ITPCsession(i) = val;
            ITPCsession(i+1) = val;
            subject = subject+1;
        end

    case 'TA2'
        sessionNames = {'R0817_20181120', 'R0817_20190625',...
            'R0898_20190723', 'R0898_20190724',...
            'R0959_20181128', 'R0959_20190703',...
            'R0983_20190722', 'R0983_20190723',...
            'R1103_20181121', 'R1103_20190710',...
            'R1187_20181119', 'R1187_20190703',...
            'R1373_20181128', 'R1373_20190708',...
            'R1452_20181119', 'R1452_20190711',...
            'R1507_20190621', 'R1507_20190627',...
            'R1547_20190729', 'R1547_20190730'}; % N=10 x 2 sessions TA2
        %   'R1535_20190708', 'R1535_20190711',...
        subjectNames = {'R0817',...
            'R0898',...
            'R0959',...
            'R0983',...
            'R1103',...
            'R1187',...
            'R1373',...
            'R1452',...
            'R1507',...
            'R1547'};
        
    otherwise
        error('expt not recognized')
end