% --- Compute Jacobian matrix for HySICS channels-----

% For the retrieval of ln(r_top), ln(r_bot), ln(tau_c)


% By Andrew J. Buggee
%%

function jacobian_ln = compute_jacobian_HySICS_logState(state_vector, measurement_estimate_ln, GN_inputs,...
    spec_response, jacobian_barPlot_flag, folder_paths)


% convert the measurement back to linear space
meas_est_linear = exp(measurement_estimate_ln);


% Define the measurement variance for the current pixel
measurement_variance = GN_inputs.measurement.variance_noLog;

% define the measurement uncertainty in linear space
measurement_uncert = GN_inputs.measurement.uncertainty(1)*100;  % percent 


% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% ----- unpack parallel for loop variables ------
% We want to avoid large broadcast variables!
wavelengths2run = GN_inputs.RT.wavelengths2run;
libRadtran_inp = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path = folder_paths.libRadtran_water_cloud_files;
mie_folder_path = folder_paths.libRadtran_mie_folder;

which_computer = GN_inputs.which_computer;




% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(GN_inputs.RT.wavelengths2run,[],"all"), max(GN_inputs.RT.wavelengths2run, [], "all")];

[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, GN_inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = GN_inputs.RT.source_file_resolution/3;   % nm



% how many wavelengths are there?
num_wl = size(GN_inputs.RT.wavelengths2run,1);

% how many state variables are there?
num_state_variables = length(state_vector);



% ---------------------------------------------------------
% ---- define the incremental change to each variable -----

% Define the fractional change the represents the partial derivative based
% on the measurement uncertainty. Define the change for each variable.
partial_diff_change = measurement_uncert.*[1/20, 1/5.7143, 1/20];
% below better for computing information content?
% change_in_state = [0.1 * r_top, 0.35 * r_bottom, 0.1 * tau_c, 0.2*wv_col_aboveCloud];
% below better for retrieval

change_in_state = partial_diff_change.* [r_top, r_bottom, tau_c];


% ----------------------------------------------------------------

% each column is a unique state vector with once variable being adjusted
state_vectors_with_change = repmat(state_vector, 1, num_state_variables) + diag(change_in_state);

% ----------------------------------------------------------

% changing variable steps through reff, tauC, and wavelength
% in for loop speak, it would be:
% for xx = 1:state_variable
%    for ww = 1:num_wl
changing_variables = [];
for xx = 1:num_state_variables

    changing_variables = [changing_variables; repmat(state_vectors_with_change(:,xx)', num_wl,1),...
        wavelengths2run];

end

% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables = [changing_variables, repmat((1:num_wl)', num_state_variables, 1)];

% how many INP files?
num_INP_files = size(changing_variables, 1);





% ----------------------------------------------------------
% create droplet profile - there are only three unique ones!
% ----------------------------------------------------------

re_with_topChange = create_droplet_profile2([changing_variables(1,1), changing_variables(1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

re_with_botChange = create_droplet_profile2([changing_variables(num_wl+1,1), changing_variables(num_wl+1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

re_with_tauChange = create_droplet_profile2([changing_variables(2*num_wl +1,1), changing_variables(2*num_wl +1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector


% -----------------------------------------------------------
% create water-cloud file - there are only three unique ones!
% -----------------------------------------------------------

wc_re_top_change = write_wc_file(re_with_topChange, changing_variables(1,3), GN_inputs.RT.z_topBottom,...
    GN_inputs.RT.lambda_forTau, GN_inputs.RT.distribution_str, GN_inputs.RT.distribution_var,...
    GN_inputs.RT.vert_homogeneous_str, GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar,...
    GN_inputs.compute_weighting_functions, which_computer, 1, 2, wc_folder_path, mie_folder_path);

wc_re_bot_change = write_wc_file(re_with_botChange, changing_variables(num_wl+1,3), GN_inputs.RT.z_topBottom,...
    GN_inputs.RT.lambda_forTau, GN_inputs.RT.distribution_str, GN_inputs.RT.distribution_var,...
    GN_inputs.RT.vert_homogeneous_str, GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar,...
    GN_inputs.compute_weighting_functions, which_computer, 2, 2, wc_folder_path, mie_folder_path);

wc_tau_change = write_wc_file(re_with_tauChange, changing_variables(2*num_wl +1,3), GN_inputs.RT.z_topBottom,...
    GN_inputs.RT.lambda_forTau, GN_inputs.RT.distribution_str, GN_inputs.RT.distribution_var,...
    GN_inputs.RT.vert_homogeneous_str, GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar,...
    GN_inputs.compute_weighting_functions, which_computer, 3, 2, wc_folder_path, mie_folder_path);

new_measurement_estimate = zeros(num_INP_files, 1);

parfor nn = 1:num_INP_files


    if nn>=1 && nn<=(num_wl)

        wc_filename = wc_re_top_change;

    elseif nn>=(num_wl+1) && nn<=(2*num_wl)

        wc_filename = wc_re_bot_change;

    elseif nn>=(2*num_wl+1) && nn<=(3*num_wl)

        wc_filename = wc_tau_change;

    end



    % define the input file name
    inputFileName = [num2str(round(mean(changing_variables(nn, end-2:end-1)))), '_','nm',...
        '_rTop_', num2str(round(changing_variables(nn,1), 4)),...
        '_rBot_', num2str(round(changing_variables(nn,2), 4)),...
        '_tauC_', num2str(round(changing_variables(nn,3), 4)), '.INP'];

    outputFileName = ['OUTPUT_',inputFileName(1:end-4)];


    % ----- Write an INP file --------
%     write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName, GN_inputs,...
%         changing_variables(nn, 4:5), wc_filename{1});
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName, GN_inputs,...
        changing_variables(nn, 4:5), wc_filename{1}, [], changing_variables(nn,3));



    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------


    % compute INP file
    runUVSPEC_ver2(libRadtran_inp, inputFileName, outputFileName, which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName, GN_inputs,...
        GN_inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(changing_variables(nn,4) - wl_perturb) &...
        source_wavelength<=(changing_variables(nn,5) + wl_perturb);

    [new_measurement_estimate(nn), ~] = reflectanceFunction_ver2(GN_inputs, ds,...
        source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



end



% Compute the change in the measurement and the jacobian matrix

% *** When transforming to the variables to log space... ***
% We transform the variables using x' = ln(x) and y' = ln(y)
% The jacobian is now: dy'/dx' = dln(y)/dln(x). This becomes
% dy'/dx' = K' = dy/dx * x/y = K * (x/y)
% ** Use equation 6.61 from Rodgers (2000). This is the same as equation C3
% from Dubovik and King (2000). This shows that the jacobian in log space
% is: dLog(F)/dLog(x) = K x/F

% Compute the change in the measurement and the jacobian matrix
change_in_measurement = reshape(new_measurement_estimate, num_wl, num_state_variables) - ...
    repmat(meas_est_linear, 1, num_state_variables);

dF_dx = change_in_measurement./repmat(change_in_state,num_wl,1);
jacobian_ln = dF_dx .* (repmat(state_vector', num_wl, 1) ./ repmat(meas_est_linear, 1, num_state_variables));

% ----- Check to see if there are any NaN values in the Jacobian Matrix -----

if any(isnan(jacobian_ln))==true

    error([newline, 'There are NaN values in the Jacobian matrix.', newline])
end



% --- Optional Plot! ---

if jacobian_barPlot_flag==true

    spectral_bands = zeros(1,length(GN_inputs.spec_response));
    for bb = 1:length(GN_inputs.spec_response)

        spectral_bands(bb) = round(median(GN_inputs.spec_response{bb}(:,1)));
    end
    [~, index_sort] = sort(spectral_bands);
    string_bands = string(spectral_bands(index_sort));


    f = figure; bar(abs(change_in_measurement(index_sort,:)))
    hold on;
    plot(sqrt(measurement_variance(index_sort)), 'k--')
    hold on
    xticklabels(string_bands);
    xlabel('Wavelength $(nm)$', 'Interpreter','latex')
    ylabel('$\triangle$ Reflectance','Interpreter','latex')
    legend('$\triangle r_{top}$','$\triangle r_{bot}$', '$\triangle \tau_{c}$','$\sigma_\lambda$',...
        'interpreter', 'latex', 'Location','best','Fontsize',20);
    grid on; grid minor
    set(f, 'Position',[0 0 1000 500])
    title('The Jacobian', 'Interpreter','latex')
    dim = [.14 0.67 .3 .3];
    str = ['$r_{top} = $',num2str(r_top),', $r_{bot} = $ ',num2str(r_bottom),', $\tau_c = $ ',num2str(tau_c)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
        'FontWeight','bold','FontSize',14, 'EdgeColor','w','Interpreter','latex');


end






end