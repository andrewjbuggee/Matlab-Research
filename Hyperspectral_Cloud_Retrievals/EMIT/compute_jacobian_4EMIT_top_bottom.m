% --- Compute Jacobian matrix for EMIT channels-----



% By Andrew J. Buggee
%%

function jacobian = compute_jacobian_4EMIT_top_bottom(emit, state_vector, measurement_estimate, inputs,...
    pixels2use, pp, jacobian_barPlot_flag)


% Define the measurement variance for the current pixel
measurement_variance = inputs.measurement.variance(:,pp);



% --- compute the Jacobian at out current estimate ---
r_top = state_vector(1);
r_bottom = state_vector(2);
tau_c = state_vector(3);


% ---------------------------------------------------------
% ---- define the incremental change to each variable -----

%change_in_state = [0.35 * r_top, 0.35 * r_bottom, 0.15 * tau_c];
%change_in_state = [0.03 * r_top, 0.25 * r_bottom, 0.04 * tau_c];        % values that just exceed measurement uncertainty for the Nov 2009 data set
change_in_state = [0.05 * r_top, 0.2 * r_bottom, 0.025 * tau_c]; 
% ----------------------------------------------------------------




% ----------------------------------------------------------

% Set up a few constants for the water cloud
H = inputs.RT.cloudDepth;                                % km - geometric thickness of cloud
n_layers = inputs.RT.cloud_layers;                          % number of layers to model within cloud

% Cloud top
z_top = inputs.RT.cloudTop_height(pp);        % km -  cloud top height

%z0 = 0.9;                                 % km - base height of cloud
z = linspace(z_top-H, z_top,n_layers);        % km - altitude above ground vector

indVar = 'altitude';                    % string that tells the code which independent variable we used

profile_type = inputs.model.profile.type; % type of water droplet profile
num_model_parameters = inputs.num_model_parameters;
dist_str = inputs.RT.drop_distribution_str;                         % droplet distribution

% -- For now, lets assume this is constant --
dist_var = linspace(inputs.RT.drop_distribution_var, inputs.RT.drop_distribution_var,...
    inputs.RT.cloud_layers);              % distribution variance
vert_homogeneous_str = inputs.RT.vert_homogeneous_str;          % This tells the function whether of not to create a multi-layered cloud
z_topBottom = [z(end), z(1)];           % km - boundaries of the altitude vector.

% Tell the code to use a pre-computed mie table for the extinction
% efficiency, or to use the value of the extinction paradox -> Qe = 2
parameterization_str = inputs.RT.parameterization_str;

% Using the same wavelength MODIS write_INP_file_4MODIS_2 uses to compute
% the cloud properties
wavelength_tau_c = emit.radiance.wavelength(inputs.bands2run(1));    % nm - Wavelength used for cloud optical depth calculation

% Lets step through each model variable and compute the derivative
jacobian = zeros(length(measurement_estimate), num_model_parameters);
change_in_measurement = zeros(length(measurement_estimate), num_model_parameters);

% ----- Let's define the 3 new state vectors -----
% each new state vector perturbs one variable only


for xx = 1:num_model_parameters
    
    % We start with the original state vector
    newState_vector = state_vector;
    
    % Then we alter just one of them
    newState_vector(xx) = state_vector(xx) + change_in_state(xx);
    
    new_r_top = newState_vector(1);
    new_r_bottom = newState_vector(2);
    new_tau_c = newState_vector(3);
    % --------------------------------------------
    % create water cloud file with droplet profile
    % --------------------------------------------
    
    new_re = create_droplet_profile2([new_r_top, new_r_bottom], z, indVar, profile_type);     % microns - effective radius vector
    
    
    loop_var = 0;

    wc_filename = write_wc_file(new_re, new_tau_c, z_topBottom, wavelength_tau_c(1,1), dist_str,...
        dist_var, vert_homogeneous_str, parameterization_str, loop_var);
    
    
    % ----- Write an INP file --------
    names.inp = write_INP_file_4EMIT_Gauss_Newton(inputs, pixels2use(pp), emit, wc_filename);
    
    % now lets write the output names
    
    names.out = writeOutputNames(names.inp);
    
    % ---- Run uvspec for the files created -----
    [new_measurement_estimate,~] = runReflectanceFunction_4EMIT_gaussNewton(names, inputs, emit.spec_response.value);
    
    change_in_measurement(:,xx) = new_measurement_estimate' - measurement_estimate;

    jacobian(:,xx) = change_in_measurement(:,xx)./change_in_state(xx);


    
end


% ----- Check to see if there are any NaN values in the Jacobian Matrix -----

if any(isnan(jacobian))==true

    error([newline, 'There are NaN values in the Jacobian matrix.', newline])
end



% --- Optional Plot! ---

if jacobian_barPlot_flag==true

    spectral_bands = zeros(1,length(inputs.bands2run));
    for bb = 1:length(inputs.bands2run)

        spectral_bands(bb) = round(median(emit.spec_response.wavelength(inputs.bands2run(bb), :)));

    end

    % turn the spectral channels into a string array
    string_bands = string(spectral_bands);

    % create a categorical array
    string_bands_cat = categorical(string_bands);
    % reorder the categorical array to fix the order of the bar chart
    string_bands_cat_reorder = reordercats(string_bands_cat, string_bands);


    f = figure; bar(string_bands_cat_reorder, abs(change_in_measurement))
    hold on;
    plot(sqrt(measurement_variance), 'k--')
    hold on

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



%-------------------------------------------------------------
% ---- SPECIAL PLOT FOR r_bot --------------------------------
%-------------------------------------------------------------

% spectral_bands = modisBands(1:7);
% [~, index_sort] = sort(spectral_bands);
% string_bands = string(round(spectral_bands(index_sort(:,1),1)));
% 
% load('jacobian_rt-10_rb-9_tau-20.mat', 'change_in_measurement')
% change_r_bot = change_in_measurement(index_sort(:,1),2);
% load('jacobian_rt-10_rb-9_tau-15.mat', 'change_in_measurement')
% change_r_bot = [change_r_bot, change_in_measurement(index_sort(:,1),2)];
% load('jacobian_rt-10_rb-9_tau-10.mat', 'change_in_measurement','measurement_variance')
% change_r_bot = [change_r_bot, change_in_measurement(index_sort(:,1),2)];
% 
% f = figure; bar([abs(change_r_bot),sqrt(measurement_variance(index_sort(:,1)))])
% hold on
% xticklabels(string_bands);
% xlabel('Wavelength $(nm)$', 'Interpreter','latex')
% ylabel('$\triangle$ Reflectance','Interpreter','latex')
% legend('$\tau_c = 20$','$\tau_c = 15$', '$\tau_c = 10$','$\sigma_\lambda$',...
%      'interpreter', 'latex', 'Location','best','Fontsize',20); 
% grid on; grid minor
% set(f, 'Position',[0 0 1000 500])
% title('$\partial F(\vec{x})/\partial r_{bot}$', 'Interpreter','latex')






end