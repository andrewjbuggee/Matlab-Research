%% Test slurm script to create Neural Network Training Data set with VOCALS-REx measurments

clear variables;

addLibRadTran_paths;

folder_paths = define_folderPaths_for_HySICS(1);

measurement_idx = 11;

mu_sample = linspace(cosd(0), cosd(65), 8);
sza = acosd(mu_sample);

% hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx);
% generate_hysics_refl_from_vocalsRex_and_era5(folder_paths, measurement_idx);
% generate_hysics_refl_from_vocalsRex_and_era5_loopGeometry(folder_paths, measurement_idx)

% *** Increase speed! now 1 par for loop instead of 2 ***
% generate_hysics_refl_from_vocalsRex_and_era5_loopGeometry_ver2(folder_paths, measurement_idx)

for nn = 1:length(sza)
    % *** 1 SZA per file so multiple can run on the supercomputer at once ***
    generate_hysics_refl_from_vocalsRex_and_era5_pick_SZA_loopGeometry_ver2(folder_paths, measurement_idx, sza(nn))
end


%% %% Test slurm script to create Neural Network Training Data set with ORACLES measurments
