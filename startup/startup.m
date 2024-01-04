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

% add generally useful functions folder
if strcmp(whatComputer, 'andrewbuggee')==true

    addpath('/Users/andrewbuggee/Documents/MATLAB/Generally_Useful_Functions')
    addpath('/Users/andrewbuggee/Documents/MATLAB/startup')

elseif strcmp(computer_name,'anbu8374')==true
% ----- MACBOOK ------

    addpath('/Users/anbu8374/Documents/MATLAB/Generally_Useful_Functions')
    addpath('/Users/anbu8374/Documents/MATLAB/startup')

elseif strcmp(computer_name,'curc')==true
% ----- CU SUPERCOMPUTER ------
    
    addpath('/projects/anbu8374/Matlab-Research/startup')

end


