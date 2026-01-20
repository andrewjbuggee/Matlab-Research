%% Test script to run the HySICS measurement function as set up on the super computer

clear variables; 

addLibRadTran_paths; 

folder_paths = define_folderPaths_for_HySICS(1); 

measurement_idx = 65; 

% hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array(folder_paths, measurement_idx);
hysics_refl_pt3_percent_in_situ_prof_and_tau_func_array_2(folder_paths, measurement_idx);

