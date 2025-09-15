
% --- Compute Jacobian matrix for HySICS channels-----


% For the retrieval of r_top, r_bot, tau_c, cwv


% By Andrew J. Buggee
%%

function jacobian = compute_jacobian_HySICS_ver2(state_vector, measurement_estimate, GN_inputs,...
    spec_response, jacobian_barPlot_flag, folder_paths)




% Define the measurement variance for the current pixel
measurement_variance = GN_inputs.measurement.variance;



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
mie_folder_path = folder_paths.libRadtran_mie_folder;
atm_folder_path = folder_paths.atm_folder_path;

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

% how many state variables are there?
num_state_variables = length(state_vector);

% ---------------------------------------------------------
% ---- define the incremental change to each variable -----

% below better for computing information content?
% change_in_state = [0.1 * r_top, 0.35 * r_bottom, 0.1 * tau_c, 0.2*wv_col_aboveCloud];
% below better for retrieval
change_in_state = [0.1 * r_top, 0.35 * r_bottom, 0.1 * tau_c, 0.1*wv_col_aboveCloud];


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
% create droplet profile - there are four unique ones!
% ----------------------------------------------------------

re_with_topChange = create_droplet_profile2([changing_variables(1,1), changing_variables(1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

re_with_botChange = create_droplet_profile2([changing_variables(num_wl+1,1), changing_variables(num_wl+1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

re_with_tauChange = create_droplet_profile2([changing_variables(2*num_wl +1,1), changing_variables(2*num_wl +1,2)],...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector

re_with_noChange = create_droplet_profile2(changing_variables(3*num_wl +1,1:2),...
    GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector


% -----------------------------------------------------------
%    create water-cloud file - there are four unique ones!
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

wc_with_no_change = write_wc_file(re_with_noChange, changing_variables(3*num_wl +1,3), GN_inputs.RT.z_topBottom,...
    GN_inputs.RT.lambda_forTau, GN_inputs.RT.distribution_str, GN_inputs.RT.distribution_var,...
    GN_inputs.RT.vert_homogeneous_str, GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar,...
    GN_inputs.compute_weighting_functions, which_computer, 4, 2, wc_folder_path, mie_folder_path);

% -----------------------------------------------------------
% create water vapor density profiles - there are only two!
% -----------------------------------------------------------
aboveCloud_waterVaporColumn_fileName_noChange = alter_aboveCloud_columnWaterVapor_profile(GN_inputs,...
    wv_col_aboveCloud, atm_folder_path);

aboveCloud_waterVaporColumn_fileName_withChange = alter_aboveCloud_columnWaterVapor_profile(GN_inputs,...
    state_vectors_with_change(end,end), atm_folder_path);


new_measurement_estimate = zeros(num_INP_files, 1);

parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files

    % initialize temporary variables
    wc_filename = {};
    waterVaporProfile_filename = [];

    % figure out which wc_filename to use
    if nn>=1 && nn<=(num_wl)

        wc_filename = wc_re_top_change;
        waterVaporProfile_filename = aboveCloud_waterVaporColumn_fileName_noChange;

        %         % define the input file name
        %         inputFileName = [num2str(mean(changing_variables(nn, 5:6))), '_','nm_rTop_', num2str(state_vectors_with_change(1,1)),...
        %             '_rBot_', num2str(r_bottom),'_tauC_', num2str(tau_c), '_CWV_', num2str(wv_col_aboveCloud),...
        %             '.INP'];


    elseif nn>=(num_wl+1) && nn<=(2*num_wl)

        wc_filename = wc_re_bot_change;
        waterVaporProfile_filename = aboveCloud_waterVaporColumn_fileName_noChange;

        %         % define the input file name
        %         inputFileName = [num2str(mean(changing_variables(nn, 5:6))), '_','nm_rTop_', num2str(r_top),...
        %             '_rBot_', num2str(state_vectors_with_change(2,2)),'_tauC_', num2str(tau_c), '_CWV_', num2str(wv_col_aboveCloud),...
        %             '.INP'];


    elseif nn>=(2*num_wl+1) && nn<=(3*num_wl)

        wc_filename = wc_tau_change;
        waterVaporProfile_filename = aboveCloud_waterVaporColumn_fileName_noChange;

        %         % define the input file name
        %         inputFileName = [num2str(mean(changing_variables(nn, 5:6))), '_','nm_rTop_', num2str(r_top),...
        %             '_rBot_', num2str(r_bottom),'_tauC_', num2str(state_vectors_with_change(3,3)), '_CWV_', num2str(wv_col_aboveCloud),...
        %             '.INP'];

    elseif nn>=(3*num_wl+1)

        wc_filename = wc_with_no_change;
        waterVaporProfile_filename = aboveCloud_waterVaporColumn_fileName_withChange;


        %         % define the input file name
        %         inputFileName = [num2str(mean(changing_variables(nn, 5:6))), '_','nm_rTop_', num2str(r_top),...
        %             '_rBot_', num2str(r_bottom),'_tauC_', num2str(tau_c), '_CWV_', num2str(state_vectors_with_change(end,end)),...
        %             '.INP'];

    end


    % define the input file name
    inputFileName = [num2str(round(mean(changing_variables(nn, end-2:end-1)))), '_','nm',...
        '_rTop_', num2str(round(changing_variables(nn,1), 4)),...
        '_rBot_', num2str(round(changing_variables(nn,2), 4)),...
        '_tauC_', num2str(round(changing_variables(nn,3), 4)),...
        '_CWV_', num2str(round(changing_variables(nn,4), 4)), '.INP'];


    outputFileName = ['OUTPUT_',inputFileName(1:end-4)];


    % ----- Write an INP file --------
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName, GN_inputs,...
        changing_variables(nn, 5:6), wc_filename{1}, [], [], waterVaporProfile_filename);



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
    idx_wl = source_wavelength>=(changing_variables(nn,5) - wl_perturb) &...
        source_wavelength<=(changing_variables(nn,6) + wl_perturb);

    [new_measurement_estimate(nn), ~] = reflectanceFunction_ver2(GN_inputs, ds,...
        source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



end



% Compute the change in the measurement and the jacobian matrix
change_in_measurement = reshape(new_measurement_estimate, num_wl, num_state_variables) - ...
    repmat(measurement_estimate, 1, num_state_variables);

jacobian = change_in_measurement./repmat(change_in_state,num_wl,1);

% ----- Check to see if there are any NaN values in the Jacobian Matrix -----

if any(isnan(jacobian))==true

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