function bp = fs_tue_basepath
% function bp = fs_tue_basepath
%
% 1. Description:
%       Base path of the Fluctuation Strength model.
% 
%   Usage: bp = amtbasepath;
%
%   `amtbasepath` returns the top level directory in which the AMT
%   files are installed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
f = mfilename('fullpath');
flast = il_strsplit(f,filesep);
L = length(flast{end})+length(flast{end-1})+1;

bp = f(1:end-L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = il_strsplit(s,d)

str={};
prev=1;
count=1;
for i=1:length(s)
    if (s(i)==d)
        if (prev<=i-1)
            str{count}=s(prev:i-1);
        end
        count=count+1;
        prev=i+1;
    end
end

str{count}=s(prev:end);