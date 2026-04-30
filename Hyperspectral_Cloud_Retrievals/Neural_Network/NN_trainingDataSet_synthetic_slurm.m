%% Test slurm script to create Neural Network Training Data set with VOCALS-REx measurments

clear variables;

addLibRadTran_paths;

folder_paths = define_folderPaths_for_HySICS(1);

measurement_idx = 1;


which_computer = whatComputer();

% ---- libRadtran folder paths ----
if strcmp(which_computer, 'andrewbuggee')

    input_file = '/Users/andrewbuggee/Documents/VS_CODE/Python-Research/lasp-CU-paper-3/training_inputs/training_inputs_jointMVN_N300000_L7.nc';

    output_dir = '/Users/andrewbuggee/Downloads/test_new_synthetic_data/';

elseif strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'curc')



end




% hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx);
% generate_hysics_refl_from_vocalsRex_and_era5(folder_paths, measurement_idx);
% generate_hysics_refl_from_vocalsRex_and_era5_loopGeometry(folder_paths, measurement_idx)

% *** Increase speed! now 1 par for loop instead of 2 ***
% generate_hysics_refl_from_vocalsRex_and_era5_loopGeometry_ver2(folder_paths, measurement_idx)

hysics_refl_from_synthetic_NN_inputs(input_file, measurement_idx, folder_paths, output_dir);
