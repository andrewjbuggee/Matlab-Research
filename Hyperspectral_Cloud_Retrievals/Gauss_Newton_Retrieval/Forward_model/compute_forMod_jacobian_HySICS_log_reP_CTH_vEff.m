
% --- Compute the forward model parameter Jacobian matrix for HySICS channels-----

% *** CURRENT UNCERTAINTIES CONSIDERED ***
% (1) Adiabatic droplet profile assumption
% (2) Cloud top height
% (3) Droplet distribution effective variance assumption


% For the retrieval of ln(r_top), ln(r_bot), ln(tau_c), and ln(acpw)


% Compute the forward model uncertainty for assuming a theoretical droplet
% profile


% By Andrew J. Buggee
%%

function jacobian_fm_ln = compute_forMod_jacobian_HySICS_log_reP_CTH_vEff(state_vector, measurement_estimate_ln, GN_inputs,...
    spec_response, jacobian_barPlot_flag, folder_paths)


% convert the measurement back to linear space
meas_est_linear = exp(measurement_estimate_ln);



% define the measurement uncertainty in linear space
measurement_uncert = GN_inputs.measurement.uncertainty(1)*100;  % percent


% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);
wv_col_aboveCloud = state_vector(4);

% ----- unpack parallel for loop variables ------
% We want to avoid large broadcast variables!
wavelengths2run = GN_inputs.RT.wavelengths2run;
libRadtran_inp = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path = folder_paths.libRadtran_water_cloud_files;
atm_folder_path = folder_paths.atm_folder_path;
mie_folder_path = folder_paths.libRadtran_mie_folder;

which_computer = GN_inputs.which_computer;




% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(wavelengths2run,[],"all"), max(wavelengths2run, [], "all")];

[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, GN_inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = GN_inputs.RT.source_file_resolution/3;   % nm



% how many wavelengths are there?
num_wl = size(GN_inputs.RT.wavelengths2run,1);


%% define the variables that are changing

% Define the length of the droplet profile, the number of layers modeled
% within the cloud
num_cloud_layers = GN_inputs.RT.n_layers;

% unpack the forward model parameters 
% (1) Adiabatic droplet profile assumption
% (2) Cloud top height assumption
% (3) Droplet distribution effective variance assumption

% The effective radius profile and the alpha parameter profile should
% change from base to top
forward_model_params = [ fliplr(GN_inputs.model.forward_model.re.mean{end}'),...
    GN_inputs.model.forward_model.cloudTopHeight.mean,...
    fliplr(GN_inputs.model.forward_model.alpha.mean{end})];


num_forward_model_params = length(forward_model_params);

% ---------------------------------------------------------
% ---- define the incremental change to each variable -----
% ---------------------------------------------------------

% Define the fractional change the represents the partial derivative based
% on the measurement uncertainty. Define the change for each variable.
partial_diff_change = measurement_uncert .*...
    [linspace(1/20, 1/5.7143, num_cloud_layers), 1/10,...
     linspace(1/10, 1/10, num_cloud_layers)];

% Compute the change in each parameter
change_in_params = partial_diff_change .* forward_model_params;


% ----------------------------------------------------------------

% each column is a unique parameter vector with one variable being adjusted
fm_params_with_change = repmat(forward_model_params', 1, num_forward_model_params) + diag(change_in_params);

% ----------------------------------------------------------

% changing variable
% in for loop speak, it would be:
% for xx = 1:state_variable
%    for ww = 1:num_wl
changing_variables = [];
for xx = 1:num_forward_model_params

    changing_variables = [changing_variables; repmat(fm_params_with_change(:,xx)', num_wl,1),...
        wavelengths2run];

end

% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables = [changing_variables, repmat((1:num_wl)', num_forward_model_params, 1)];



% -----------------------------------------------------------
% ---------- create water vapor density profiles ------------
% -----------------------------------------------------------
% ** create file with original cloud top height **
aboveCloud_waterVaporColumn_fileName_original_cloudTopHeight = alter_aboveCloud_columnWaterVapor_profile(GN_inputs,...
    wv_col_aboveCloud, atm_folder_path);

% ** create file with new cloud top height **
% Define the cloud top value
GN_inputs.RT.z_topBottom = [changing_variables(num_forward_model_params * num_wl, num_cloud_layers+1), ...
    changing_variables(num_forward_model_params * num_wl, num_cloud_layers+1) - GN_inputs.RT.cloud_depth];

aboveCloud_waterVaporColumn_fileName_new_cloudTopHeight = alter_aboveCloud_columnWaterVapor_profile(GN_inputs,...
    wv_col_aboveCloud, atm_folder_path);
% -----------------------------------------------------------
% -----------------------------------------------------------


% how many INP files?
num_INP_files = size(changing_variables, 1);


% create a cell array for ACPW filenames
acpw_filenames = cell(num_INP_files, 1);
% Define the ACPW filenames
acpw_filenames(1:(num_cloud_layers * num_wl)) =  {aboveCloud_waterVaporColumn_fileName_original_cloudTopHeight};
acpw_filenames((num_cloud_layers * num_wl)+1 : end) =  {aboveCloud_waterVaporColumn_fileName_new_cloudTopHeight};




wc_filenames = cell(num_forward_model_params, 1);

for xx = 1:num_forward_model_params

        % --------------------------------------
        % ** Set the effective radius profile **

        % create droplet profile
        re_profile = changing_variables(xx, 1:num_cloud_layers);

        % for the function 'write_wc_file', the re vector must be arranged so
        % that the first entry is at cloud bottom, and the final value is at
        % cloud top. Flip the above vector
        re_profile = fliplr(re_profile);
        % --------------------------------------


        % --------------------------------------
        % ****** Set the cloud top height ******
        GN_inputs.RT.z_topBottom = [changing_variables(xx, num_cloud_layers+1), ...
            changing_variables(xx * num_wl, num_cloud_layers+1) - GN_inputs.RT.cloud_depth];
        % --------------------------------------
    
        
        % -----------------------------------------
        % **** Set the alpha parameter profile ****
        % for the function 'write_wc_file', the alpha vector must be arranged so
        % that the first entry is at cloud bottom, and the final value is at
        % cloud top.
        alpha_profile = changing_variables(xx, (num_cloud_layers+2):((2*num_cloud_layers+1)));
        % -----------------------------------------


    % -----------------------------------------------------------
    % ---------------- create water-cloud file ------------------
    % -----------------------------------------------------------

    wc_filenames{xx} = write_wc_file(re_profile, tau_c, GN_inputs.RT.z_topBottom,...
        GN_inputs.RT.lambda_forTau, GN_inputs.RT.distribution_str, alpha_profile,...
        GN_inputs.RT.vert_homogeneous_str, GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar,...
        GN_inputs.compute_weighting_functions, which_computer, xx, 2, wc_folder_path, mie_folder_path);

end



% repeat the wc_filenames array for easy access in the for loop below
wc_filenames = reshape( repmat(wc_filenames', num_wl, 1), num_INP_files, 1);



new_measurement_estimate = zeros(num_INP_files, 1);




parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files

    % grab the paramter number
    param_str = extractBetween(wc_filenames{nn}, '_nn', '.DAT');

    % define the input file name
    inputFileName = [num2str(round(mean(changing_variables(nn, end-2:end-1)))), '_','nm',...
        'fm_jacobian_parameter_', param_str{1}, '.INP'];


    outputFileName = ['OUTPUT_',inputFileName(1:end-4)];


    % ----- Write an INP file --------
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName, GN_inputs,...
        changing_variables(nn, end-2:end-1), wc_filenames{nn}{1}, [], tau_c,...
        acpw_filenames{nn});



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
    idx_wl = source_wavelength>=(changing_variables(nn,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables(nn, end-1) + wl_perturb);

    [new_measurement_estimate(nn), ~] = reflectanceFunction_ver2(GN_inputs, ds,...
        source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



end



% *** When transforming to the variables to log space... ***
% We transform the variables using x' = ln(x) and y' = ln(y)
% The jacobian is now: dy'/dx' = dln(y)/dln(x). This becomes
% dy'/dx' = K' = dy/dx * x/y = K * (x/y)
% ** Use equation 6.61 from Rodgers (2000). This is the same as equation C3
% from Dubovik and King (2000). This shows that the jacobian in log space
% is: dLog(F)/dLog(x) = K x/F

% Compute the change in the measurement and the jacobian matrix
change_in_measurement = reshape(new_measurement_estimate, num_wl, num_forward_model_params) - ...
    repmat(meas_est_linear, 1, num_forward_model_params);

dF_dx = change_in_measurement./repmat(change_in_params, num_wl,1);
jacobian_fm_ln = dF_dx .* (repmat(forward_model_params, num_wl, 1) ./ repmat(meas_est_linear, 1, num_forward_model_params));

% ----- Check to see if there are any NaN values in the Jacobian Matrix -----

if any(isnan(jacobian_fm_ln), 'all')==true

    error([newline, 'There are NaN values in the Jacobian matrix.', newline])
end



% --- Optional Plot! ---

if jacobian_barPlot_flag==true

    % Define the measurement variance for the current pixel
    measurement_variance = GN_inputs.measurement.variance_noLog;

    spectral_bands = zeros(1,length(GN_inputs.spec_response));
    for bb = 1:length(GN_inputs.spec_response)

        spectral_bands(bb) = round(median(GN_inputs.spec_response{bb}(:,1)));
    end
    [~, index_sort] = sort(spectral_bands);
    string_bands = string(spectral_bands(index_sort));


    f = figure; bar(abs(change_in_measurement(index_sort,:)))
    hold on;
    plot(sqrt(measurement_variance_ln(index_sort)), 'k--')
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