%% Custom Startup File which alters MATLAB's default settings

% By Andrew J. Buggee
% --------------------

%% ----- Setting New Default Plot Values -----

% --- Text Settings ---
set(groot,'DefaultTextColor','w');

% --- Axes Settings ---
set(groot,'DefaultAxesColor','k')
set(groot,'DefaultAxesXColor','w')
set(groot,'DefaultAxesYColor','w')
set(groot,'DefaultAxesZColor','w')
set(groot,'DefaultAxesGridColor','w')
set(groot,'DefaultAxesMinorGridColor','w')
set(groot,'DefaultAxesFontSize',18);
set(groot,'DefaultAxesFontWeight','bold');
set(groot,'DefaultAxesTitleFontSizeMultiplier',1.8);
set(groot,'DefaultAxesLabelFontSizeMultiplier',1.6);

% --- Line Settings ---
set(groot,'DefaultLineLineWidth',4);
set(groot,'DefaultLineMarkerSize',15)

% --- Figure Settings ---
set(groot,'DefaultFigureColor','k')

%% ---- Add these folders to your path -----

% First, find out which computer is being used
[status,username] = system('whoami');

username = username(1:end-1);

% determine if we're on the CURC supercomputer or my LASP mac desktop
if strcmp(username, 'anbu8374')==true && strcmp(matlabroot, '/Applications/MATLAB_R2022b.app')==true

    % Then we're using my Mac desktop at LASP
    % Keep this username!

elseif strcmp(username, 'anbu8374')==true && strcmp(matlabroot, '/scratch/local/MATLAB')==true

    % Then were on the super computer! Change the username to reflect this
    username = 'curc';

end



if status ~= 0
    error(['Status of command returned value of ',num2str(status)])
end




% add generally useful functions folder and the startup folder
if strcmp(username, 'andrewbuggee')==true

    addpath('/Users/andrewbuggee/Documents/MATLAB/Generally_Useful_Functions')
    addpath('/Users/andrewbuggee/Documents/MATLAB/startup')

elseif strcmp(username,'anbu8374')==true
% ----- MACBOOK ------

    addpath('/Users/anbu8374/Documents/MATLAB/Generally_Useful_Functions')
    addpath('/Users/anbu8374/Documents/MATLAB/startup')

elseif strcmp(username,'curc')==true
% ----- CU SUPERCOMPUTER ------
    
    addpath('/projects/anbu8374/Matlab-Research/startup')

end


