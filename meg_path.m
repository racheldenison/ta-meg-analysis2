function meg_path(option)

% function meg_path(option)
%
% option = 'add' or 'remove'. default is 'add'.
% adds or removes path to toolbox
% run from the top-level toolbox directory
%
% Rachel Denison
% March 2020

if nargin==0
    option = 'add';
end

currdir = pwd;
d = dir;

% find all subdirectories (exclude hidden directories)
subdirs = {};
for i = 1:numel(d)
    if d(i).isdir && ~strcmp(d(i).name(1),'.')
        subdirs = [subdirs, {d(i).name}];
    end
end
n = numel(subdirs);

p = '';
for i = 1:n
    subdir = subdirs{i};
    if ~isdir(subdir)
        error('%s is not a directory', subdir)
    end
    p = sprintf('%s%s/%s:', p, currdir, subdir);
end

switch option
    case 'add'
        addpath(p)
    case 'remove'
        rmpath(p)
    otherwise
        error('path option not recognized')
end

