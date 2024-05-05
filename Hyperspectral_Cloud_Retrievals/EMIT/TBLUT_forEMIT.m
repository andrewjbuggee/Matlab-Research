% Compute Cloud Optical Properties from EMIT data using the two-wavelength
% look-up-table method

% By Andrew John Buggee

%%

function [opt_prop] = TBLUT_forEMIT(emit_)

% Check to make sure we have the correct inputs

%% what computer are we using?

% a template file has been set up to be edited and saved as a new file
% determine which computer is being used
userName = whatComputer;

if strcmp(userName,'anbu8374')

    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];

elseif strcmp(userName,'andrewbuggee')

    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];

else
    error('I dont recognize this computer user name')
end


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_modis_inputs(folderName, L1B_1km_fileName);

% *** Check Inputs ***


end