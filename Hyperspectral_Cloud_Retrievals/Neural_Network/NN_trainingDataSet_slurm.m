%% Test slurm script to create Neural Network Training Data set with VOCALS-REx measurments

clear variables; 

addLibRadTran_paths; 

folder_paths = define_folderPaths_for_HySICS(1); 

measurement_idx = 65; 

% hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx);
generate_hysics_refl_from_vocalsRex_and_era5(folder_paths, measurement_idx);



%% %% Test slurm script to create Neural Network Training Data set with ORACLES measurments
