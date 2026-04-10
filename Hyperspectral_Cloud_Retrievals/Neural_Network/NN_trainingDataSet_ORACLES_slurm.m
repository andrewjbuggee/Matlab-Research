%% Local test script for generating a Neural Network training data set
%  from ORACLES in-situ ensemble profiles.
%
%  Mirrors NN_trainingDataSet_slurm.m (VOCALS-REx version).
%  On the CU Alpine supercomputer this script is driven by a SLURM shell
%  script that sets measurement_idx and output_dir via environment variables
%  or hard-coded arguments; here those values are set manually for local
%  testing.

clear variables;

addLibRadTran_paths;

folder_paths = define_folderPaths_for_HySICS(1);

% -------------------------------------------------------------------------
% --- Set the ORACLES ensemble profile index to process -------------------
% -------------------------------------------------------------------------
% Run 'length(ensemble_profiles)' after loading the ensemble file to find
% the valid index range.  The saved ensemble_profiles_with_precip_from_33_files
% file (13-Mar-2026) contains 243 profiles.

measurement_idx = 106;

% -------------------------------------------------------------------------
% --- Define the solar zenith angle grid ----------------------------------
% -------------------------------------------------------------------------
% Sample cos(sza) linearly so the sampling is uniform in optical path-length
% space.  8 values from sza=0 to sza=65 degrees.

mu_sample = linspace(cosd(0), cosd(65), 8);
sza = acosd(mu_sample);     % degrees

% -------------------------------------------------------------------------
% --- Output directory (adjust for local testing) -------------------------
% -------------------------------------------------------------------------
% On the supercomputer, output_dir is passed by the SLURM wrapper.
% Uncomment and edit the line below for local testing:

which_computer = whatComputer();

if strcmp(which_computer, 'anbu8374')

    output_dir = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/Neural_Network/Training_data_set/', ...
        'ORACLES_test/'];

elseif strcmp(which_computer, 'andrewbuggee')

    output_dir = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/Neural_Network/Training_data_set/', ...
        'ORACLES_test/'];

elseif strcmp(which_computer, 'curc')

    % On Alpine the SLURM script typically exports output_dir; define a
    % fallback here in case it is not already set.
    if ~exist('output_dir', 'var')
        output_dir = ['/scratch/alpine/anbu8374/neural_network_training_data/', ...
            'ORACLES/'];
    end

end

% -------------------------------------------------------------------------
% --- Loop over SZA values (one call per file for parallel SLURM jobs) ----
% -------------------------------------------------------------------------

for nn = 1:length(sza)
    hysics_refl_from_oracles_and_era5_SZA_loopGeometry(folder_paths, measurement_idx, sza(nn), output_dir)
end
