% --- Compute Jacobian matrix for HySICS channels-----



% By Andrew J. Buggee
%%

function jacobian = compute_jacobian_HySICS(hysics, state_vector,measurement_estimate,GN_inputs,...
    jacobian_barPlot_flag, folder_paths)


% Define the measurement variance for the current pixel
measurement_variance = GN_inputs.measurement.variance;

% --- Define the filename to save all calculations ---
saveCalculations_fileName = folder_paths.HySICS_retrievals;

% --- Define the INP Folder location ---
INP_folderName = folder_paths.libRadtran_inp;

% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% how many wavelengths are there?
num_wl = size(GN_inputs.RT.wavelengths2run,1);

% how many state variables are there?
num_state_variables = length(state_vector);

% ---------------------------------------------------------
% ---- define the incremental change to each variable -----

change_in_state = diag([0.1 * r_top, 0.35 * r_bottom, 0.1 * tau_c]); 
% ----------------------------------------------------------------

% each column is a unique state vector with once variable being adjusted
state_vectors_with_change = repmat(state_vector, 1, num_state_variables) + change_in_state;

% ----------------------------------------------------------

% changing variable steps through reff, tauC, and wavelength
% in for loop speak, it would be:
% for xx = 1:stat_variable
%    for ww = 1:num_wl
changing_variables = [];
for xx = 1:num_state_variables

changing_variables = [changing_variables; repmat(state_vectors_with_change(:,xx)', num_wl,1),...
    GN_inputs.RT.wavelengths2run];

end

% Add a final column that includes the index for the spectral response
% function. These always increase chronologically
changing_variables = [changing_variables, repmat((1:num_wl)', num_state_variables, 1)];

% how many INP files?
num_INP_files = size(changing_variables, 1);


% Lets step through each model variable and compute the derivative
jacobian = zeros(length(measurement_estimate),num_model_parameters);
change_in_measurement = zeros(length(measurement_estimate),num_model_parameters);


% ----- Let's define the 3 new state vectors -----
% each new state vector perturbs one variable only
%perturbed_state_vector = repmat(state_vector, length(state_vector), 1) + change_in_state;

% make an empty strucutre for the file names
%names = struct([]);

for nn = 1:num_INP_files
    

    % --------------------------------------------
    % create water cloud file with droplet profile
    % --------------------------------------------
    
    new_re = create_droplet_profile2([changing_variables(nn,1), changing_variables(nn,2)],...
        GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector
    
    
    loop_var = 0;

    wc_filename = write_wc_file(new_re, new_tau_c, z_topBottom, wavelength_tau_c(1,1), dist_str,...
        dist_var, vert_homogeneous_str, parameterization_str, loop_var);
    
    
    % ----- Write an INP file --------
    names.inp = write_INP_file_4MODIS_Gauss_Newton(GN_inputs, modisInputs, pixel_row, pixel_col, modis, wc_filename);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
    
    % ---- Run uvspec for the files created -----
    [new_measurement_estimate,~] = runReflectanceFunction_4modis_gaussNewton(names,INP_folderName,saveCalculations_fileName, GN_inputs.spec_response);
    
    change_in_measurement(:,xx) = new_measurement_estimate' - measurement_estimate;

    jacobian(:,xx) = change_in_measurement(:,xx)./change_in_state(xx);


    
end


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