%% Estimate the above cloud precipitable water amount using Simualted HySICS measurements


% By Andrew John Buggee

%%

function [acpw_retrieval] = ACPW_retrieval_for_HySICS(simulated_measurements, tblut_retrieval, folder_paths, print_status_updates, print_libRadtran_err)



%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_acpw = create_HySICS_inputs_ACPW(simulated_measurements.inputs, tblut_retrieval, print_libRadtran_err);


%% Find the measurements closest to the bands to run

[~, idx_1] = min(abs(simulated_measurements.inputs.bands2run - inputs_acpw.bands2run(1)));

[~, idx_2] = min(abs(simulated_measurements.inputs.bands2run - inputs_acpw.bands2run(2)));

[~, idx_3] = min(abs(simulated_measurements.inputs.bands2run - inputs_acpw.bands2run(3)));


% error if the values found are at least 15nm from the intended wavlengths
if abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_1,:)) - 872)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 872 nm', newline])

elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_2,:)) - 900)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 900 nm', newline])


elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_3,:)) - 1127)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 900 nm', newline])

else

    % Then we set the bands to run to be to ones found to be closest to the
    % desired bands out of the measurement bands provided
    inputs_acpw.bands2run_from_set_of_measurements = [idx_1, idx_2, idx_3];
    inputs_acpw.bands2plot = inputs_acpw.bands2run;

    % ---- Define the wavelengths ----
    inputs_acpw.RT.wavelengths2run = simulated_measurements.inputs.RT.wavelengths2run(inputs_acpw.bands2run_from_set_of_measurements,:);



end

%% Define the reflectance measurements used for the retrieval


% The window channel where water vapor absorption is negligible
R_window_872 = simulated_measurements.Refl_model(8);        % 1/sr - HySICS channel centered around 881 nm

% The weakly absorbing water vapor absorption channel
R_waterVap_900 = simulated_measurements.Refl_model(10);     % 1/sr - HySICS channel centered around 900 nm

% The strongly absorbing water vapor absorption channel
R_waterVap_1127 = simulated_measurements.Refl_model(21);     % 1/sr - HySICS channel centered around 1127 nm

R_measured = [R_window_872; R_waterVap_900; R_waterVap_1127];




end