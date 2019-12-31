function FS_TUe_start
% function FS_TUe_start
%
% 1. Description: 
%       - Adds to path the required folders to run the fluctuation strength
%         model.
%
% 2. Stand-alone example:
%       FS_TUe_start;
%
% 4. Additional info:
%   Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2018
%       (ale.a.osses@gmail.com), HearingTechnology, UGent, Belgium 2018-2020
% Created on    : 31/12/2019
% Last edited on: 31/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirmain  = [cd filesep];
bSuccess = isdir([dirmain 'Utility' filesep]);

if bSuccess == 0
    error('Please set the current MATLAB directory to the location of the FluctuationStrength_TUe model and re-run this script');
end

subdirs = il_get_AM_AddOns_paths(dirmain);
for i = 1:length(subdirs)
  addpath(subdirs{i})
end
dir = [dirmain 'Utility' filesep]; addpath(dir);

disp(['EOF: ' mfilename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subdirs = il_get_AM_AddOns_paths(dir)

subdirs{1} = [dir 'demos'       filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%