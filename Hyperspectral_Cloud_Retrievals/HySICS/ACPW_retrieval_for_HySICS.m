%% Estimate the above cloud precipitable water amount using Simualted HySICS measurements


% By Andrew John Buggee

%%

function [acpw_retrieval] = ACPW_retrieval_for_HySICS(simulated_measurements, folder_paths, print_status_updates, print_libRadtran_err)



%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_acpw = create_HySICS_inputs_ACPW(folder_paths, simulated_measurements.inputs, print_libRadtran_err);


%% Define the reflectance measurements used for the retrieval


% The window channel where water vapor absorption is negligible
R_window_872 = simulated_measurements.Refl_model(8);        % 1/sr - HySICS channel centered around 881 nm

% The weakly absorbing water vapor absorption channel
R_waterVap_900 = simulated_measurements.Refl_model(10);     % 1/sr - HySICS channel centered around 900 nm

% The strongly absorbing water vapor absorption channel
R_waterVap_1127 = simulated_measurements.Refl_model(21);     % 1/sr - HySICS channel centered around 1127 nm

R_measured = [R_window_872; R_waterVap_900; R_waterVap_1127];




end