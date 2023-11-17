%% ----- Read MODIS Spectral Bins -----

% This code will provide use with the modis spectral bins. They are
% organized by number. The ouput will be the upper and lower bounds, as
% well as the center wavelength. All output values are specified in
% nanometers

% inputs:
%   - n: the spectral bin number you wish to have. Please provide a vector
%   with integers, where n can be any number between 1 and 36.

% By Andrew J. Buggee
%%

function [mod_bands] = readMODISbands(n)
    
% -------------------------------------------------
% Check inputs 

    if isnumeric(n)
        
    else
        error('input must be a numeric entry')
        
    end
    
    if length(n)>36
        error('MODIS has 36 spectral bands. You requested more than 36. You may pull any of of them')
    end
    
% --------------------------------------------------

% check which computer you're using

comp = whatComputer;

% pull desired bins

if strcmp(comp,'andrewbuggee')==true
    
    folder_bands = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/MODIS_Cloud_Retrieval/MODIS_data/';
    
elseif strcmp(comp,'anbu8374')==true
    
    folder_bands = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/';
    
end

    file_bands = 'modis_bands.txt';
    fileID = fopen([folder_bands,file_bands],'r');
    T = fscanf(fileID,'%d %d %d\n');
    
    bandNumber_index = 1:3:length(T);
    bandNumber = T(bandNumber_index);
    
    % retrieve the lower bound
    bandLow_index = 2:3:length(T);
    bandLow = T(bandLow_index);
    
    % retrieve the upper bound
    bandUpper_index = 3:3:length(T);
    bandHigh = T(bandUpper_index);
    
    
    mod_bands.number = n;
    mod_bands.lowerBound = bandLow(n);
    mod_bands.upperBound = bandHigh(n);
    mod_bands.center = (mod_bands.upperBound - mod_bands.lowerBound)/2 + mod_bands.lowerBound;
    
end
