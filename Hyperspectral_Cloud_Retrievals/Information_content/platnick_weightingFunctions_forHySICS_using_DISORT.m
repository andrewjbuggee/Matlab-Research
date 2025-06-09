%% Compute weighting functions using DISORT with a high number of streams

% By Andrew John Buggee

clear variables


%% Which computer are you using?


clear inputs

% Determine which computer this is being run on
inputs.which_computer = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(inputs.which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/weighting_functions/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/weighting_functions/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = '/scratch/alpine/anbu8374/weighting_functions/';



    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = '/scratch/alpine/anbu8374/weighting_functions/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_inp, 'dir')

    mkdir(inputs.folderpath_inp)

end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_reflectance, 'dir')

    mkdir(inputs.folderpath_reflectance)

end





%%  Delete old files?
% First, delete files in the HySICS folder
delete([inputs.folderpath_inp, '*.INP'])
delete([inputs.folderpath_inp, '*.OUT'])




%% set up the inputs to create an INP file for DISORT!

[inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, false);


%% Write the INP files


% num wavelengths
num_wl = size(inputs.RT.wavelengths2run,1);


tic

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------









elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------



    % length of each independent variable
    % length of each independent variable
    % num_rTop = length(inputs.RT.r_top);
    % num_rBot = length(inputs.RT.r_bot);
    % num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);
    num_tau_layers = inputs.RT.n_layers;

    if inputs.compute_weighting_functions==true
        num_INP_files = num_wl*num_tau_layers;
    else
        num_INP_files = num_wl;
    end


    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);



    % need an re and tau_c file for each layer, because we slowly build the
    % cloud but adding layer after layer
    idx_unique_indVars = 1:num_tau_layers;


    % create a droplet profile
    re = create_droplet_profile2([inputs.RT.r_top, inputs.RT.r_bot], inputs.RT.z,...
        inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector

    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------
    % re must be defined from cloud bottom to cloud top
    % z_topBottom must be defined as the cloud top height first, then
    % cloud bottom
    [wc_filename, lwc, ext_bulk_coeff_per_LWC] = write_wc_file(re, inputs.RT.tau_c,...
        inputs.RT.z_topBottom, inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
        inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
        inputs.RT.indVar, inputs.compute_weighting_functions, inputs.which_computer, 1);


    % Compute the optical depth of each layer
    % lwc and ext_bulk_coeff_per_LWC are reported from cloud bottom to
    % cloud top
    % therefore, the optical depth of each layer starts at cloud bottom
    % inputs.RT.tau_layers = (lwc.*ext_bulk_coeff_per_LWC.*(inputs.RT.z_edges(2) - inputs.RT.z_edges(1)))';  % the optical depth of each layer, starting from cloud top
    inputs.RT.tau_layers = (lwc.*ext_bulk_coeff_per_LWC.* diff(inputs.RT.z_edges'));  % the optical depth of each layer, starting from cloud top




    % Create the changing variables matrix that defines how each INP file
    % is unique
    if inputs.RT.monochromatic_calc==true

        % Variables by column:
        % (1) - monochromatic wavelength to run
        % (2) - total optical depth (sum of all layers)
        % (3) - index indicating which spectral response function is needed
        changing_variables = [reshape(repmat(inputs.RT.wavelengths2run(:,1)', num_tau_layers,1), [],1)...
            repmat(flipud(cumsum(flipud(inputs.RT.tau_layers))), num_wl,1)];

    else

        changing_variables = [repmat(inputs.RT.wavelengths2run, num_tau_layers,1),...
            repmat(flipud(cumsum(flipud(inputs.RT.tau_layers'))), num_wl,1)];

    end


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, reshape(repmat((1:num_wl), num_tau_layers, 1), [],1)];


    % the wc_filenames should be the same for different wavelengths
    wc_filename = repmat(wc_filename, num_wl, 1);


    % Now write all the INP files
    parfor nn = 1:num_INP_files
        % for nn = 1:num_INP_files




        if inputs.RT.monochromatic_calc==true

            % set the wavelengths for each file
            wavelengths = changing_variables(nn,1);

            % ------------------------------------------------
            % ---- Define the input and output filenames! ----
            % ------------------------------------------------
            % input_names need a unique identifier. Let's give them the nn value so
            % they can be traced, and are writen over in memory


            inputFileName{nn} = ['weightingFunction_',num2str(mean(wavelengths)), '_','nm_rTop_', num2str(inputs.RT.r_top),...
                '_rBot_', num2str(inputs.RT.r_bot),'_tauC_', num2str(round(changing_variables(nn,2),4)),...
                '_layers1-', num2str(mod(num_INP_files - (nn), num_tau_layers)+1), '.INP'];



            outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        else

            % set the wavelengths for each file
            wavelengths = changing_variables(nn, 1:2);

            % ------------------------------------------------
            % ---- Define the input and output filenames! ----
            % ------------------------------------------------
            % input_names need a unique identifier. Let's give them the nn value so
            % they can be traced, and are writen over in memory


            inputFileName{nn} = ['monteCarlo_',num2str(mean(wavelengths)), '_','nm_rTop_', num2str(inputs.RT.r_top),...
                '_rBot_', num2str(inputs.RT.r_bot),'_tauC_', num2str(round(changing_variables(nn,3),4)),...
                '_layers1-', num2str(num_INP_files - (nn-1)), '.INP'];



            outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        end




        % ------------------ Write the INP File --------------------
        write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn});


    end



end

toc



%% Compute Reflectance

% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
if isscalar(inputs.RT.sensor_altitude)
    Refl_model = zeros(num_INP_files, 1);
elseif ischar(inputs.RT.sensor_altitude)
    Refl_model = zeros(num_INP_files, 1);
elseif length(inputs.RT.sensor_altitude)>1
    Refl_model = zeros(num_INP_files, length(inputs.RT.sensor_altitude));
end


parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


    disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------


    % compute INP file
    [inputSettings] = runUVSPEC(inputs.folderpath_inp, inputFileName{nn}, outputFileName{nn},...
        inputs.which_computer);

    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC(inputs.folderpath_inp, outputFileName{nn},inputSettings(2,:),...
        inputs.RT.compute_reflectivity_uvSpec);

    if inputs.RT.compute_reflectivity_uvSpec==true

        % Store the reflectance
        % Refl_model(nn, :) = ds.reflectivity.value;       % reflectance is in units of 1/sr

        % Store the azimuthally averaged reflectance
        % Azimuthal average at this particular zenith angle
        Refl_model(nn, :) = ds.reflectivity.az_avg;       % reflectance is in units of 1/sr

    else

        % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
        [Refl_model(nn, :), ~] = reflectanceFunction(inputSettings(2,:), ds,...
            spec_response_value(changing_variables(nn,end),:));

    end



end






toc




%% Compute the weighting functions using Platnick (2000)
% Equation 4 defines the weighting function as the normalized derivative of
% reflectance with respect to optical depth.

if inputs.RT.monochromatic_calc==true

    % reshape Refl_model
    Refl_model = reshape(Refl_model, inputs.RT.n_layers, size(inputs.RT.wavelengths2run,1));

    % reshape the optical depth of each file
    tau = reshape(changing_variables(:,2), inputs.RT.n_layers, size(inputs.RT.wavelengths2run,1));

    % compute the derivative of reflectivity as a function of optical depth
    % normalize by the reflectance over the full cloud optical thickness
    w = diff(flipud(Refl_model))./diff(flipud(tau)) ./ repmat(Refl_model(1,:), inputs.RT.n_layers-1, 1);

else

    % compute the derivative of reflectivity as a function of optical depth
    w = diff(flipud(Refl_model))./diff(flipud(changing_variables(:,3)));

end




%% plot the weighting function

figure;



if inputs.RT.monochromatic_calc==true

    tau_midPoint = flipud(tau);
    tau_midPoint = tau_midPoint(1:end-1,:) + diff(tau_midPoint, 1, 1);
    % plot(w, diff(flipud(changing_variables(:,2)))/2 + flipud(changing_variables(2:end,2)))
    % Plot the tau mid-point of each layer
    plot(w, tau_midPoint)

    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    %f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);

else

    tau = changing_variables(2:end,3);

    % plot(w, diff(flipud(changing_variables(:,2)))/2 + flipud(changing_variables(2:end,2)))
    plot(w, flipud(tau))

    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);

end




% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$w(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
title(['Weighting Function at ', num2str(changing_variables(1,1)), ' nm'],'Interpreter','latex')




set(gcf, 'Position',[0 0 1400 800])



% Plot the spline fit
% hold on
% tau_2plot = (0.5:0.1:inputs.RT.tau_c);
% plot(f(tau_2plot),tau_2plot, 'Color', mySavedColors(1, 'fixed'))


% Create Legend
legend(string(inputs.RT.wavelengths2run(:,1))','Interpreter','latex','Location','northwest','FontSize',22)



%% Let's renormalize the weighting functions so that they integrate to 1

for ww = 1:size(inputs.RT.wavelengths2run, 1)

    a = 1/trapz(tau_midPoint(:,ww), w(:,ww));

    w(:,ww) = w(:,ww).*a;

end


%%
% ----------------------------------------------
% ---------- SAVE REFLECTANCE OUTPUT! ----------
% ----------------------------------------------

% Save the version without an measurement uncertainty. Then we can add
% uncertainty and save the new file

if strcmp(inputs.which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    inputs.folderpath_2save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/weighting_functions/'];



elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/weighting_functions/'];


elseif strcmp(inputs.which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    warning([newline, 'No folder to store things in!', newline])



end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_2save, 'dir')

    mkdir(inputs.folderpath_2save)

end



rev = 1;

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true


    filename = [inputs.folderpath_2save,'disort_HySICS_reflectance_for_weightingFunctions_',...
        'inhomogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];


else

    filename = [inputs.folderpath_2save,'disort_HySICS_reflectance_for_weightingFunctions_',...
        'homogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

end


while isfile(filename)
    rev = rev+1;
    if rev<10
        filename = [filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        filename = [filename(1:end-6), num2str(rev),'.mat'];
    end
end


save(filename, "Refl_model","inputs", "spec_response", "changing_variables", "w", "tau_midPoint");



%% Construct the Platnick (2000) covriance matrix


% First, step through each weighting function and compute the estimated
% single band effective radius retrieval (Platnick (2000), eq. 3)
% create a droplet profile
re = create_droplet_profile2([inputs.RT.r_top, inputs.RT.r_bot], inputs.RT.z,...
    inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector

re_midPoint = flipud(re(1:end-1) + diff(re)/2);

re_est = zeros(1, size(w,2));

for ww = 1:size(w,2)

    re_est(ww) = trapz(tau_midPoint(:,ww), re_midPoint .* w(:,ww));

end


% The elements of this symmetric matrix are given by the inner product of
% the weightings (Platnick (2000) pg. 22933)

C = zeros(size(inputs.RT.wavelengths2run,1), size(inputs.RT.wavelengths2run,1));

% this matrix is symmetric, meaning C(1,3) = C(3,1).

% Normalize the covariance matrix by the two effective radii estimates for
% the ith and jth element

for ii = 1:size(inputs.RT.wavelengths2run,1)
    for jj = 1:size(inputs.RT.wavelengths2run,1)

        C(ii,jj) = trapz(tau_midPoint(:, ii), w(:,ii) .* w(:,jj))/(re_est(ii)*re_est(jj));

    end
end


% compute the eigenvalues of our scaled covariance matrix
eig_val = eig(C);


% According to Platnick (2000): "Then for all weightings to contribute unique 
% information,the minimum eigenvalue of the scaled covariance matrix must be
% greater than about e^2/(N r_m^2), where N is the number of measurements,
% r_m is the mean value of the unknown re(tau) profile, and e is the relative
% uncertainty of the retrieval and the weighting function (forward model)
% [Twomey, 1977]

% let's assume a relative uncertainty of 5%
relative_uncert = 0.05;

twomey_fraction = relative_uncert^2 / (size(w,2) * mean(re)^2);


if min(eig_val)>twomey_fraction

    disp([newline, 'All wavelengths contribute unique information!', newline])

else

    % How many unique pieces of information do these weightings contribute?
    
end



%% Plot the estimated effective radii versus wavelength and versus optical depth

% First, plot each estimated effective radii versus wavelength
figure;

subplot(1,2,1)

plot(inputs.RT.wavelengths2run(:,1), re_est, '.', 'markersize', 15)
xlabel('Wavelength (nm)','FontWeight','bold','Interpreter','latex', 'Fontsize', 35)
ylabel('$r_{e}^{*}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);
grid on; grid minor


% overlay the reflectance over all wavelengths for a cloud scene

if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    folderpath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_ALL_BANDS.mat';


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the Simulated HySICS data folder path

    folderpath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';

    % filename
    filename = 'simulated_measurement_HySICS_reflectance_inhomogeneous_droplet_profile_sim-ran-on-02-Jun-2025_rev1.mat';


end

sm = load([folderpath, filename]);

hold on
yyaxis right
plot(sm.inputs.RT.wavelengths2run(:,1), sm.Refl_model, 'linewidth', 1.5)

ylabel('Reflectance $1/sr$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);


% On the other side, plot the true droplet profile
subplot(1,2,2)
plot(re, tau(:,1))
% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$r_e$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% let's also plot every unique re retrieval within 0.1 microns
% round the estimated retrieved effective radii to the nearest 0.1 microns
% and find all unique values
re_2Plot = unique(round(re_est, 1));

% find the indexes for each re_2Plot that align with their associated tau
% level in the droplet profile
idx = zeros(length(re_2Plot), 1);
for nn = 1:length(re_2Plot)

    [~, idx(nn)] = min(abs(re - re_2Plot(nn)));

    hold on
    yline(tau(idx(nn), 1), 'LineWidth', 1, 'Color', 'k')
    
end



% set figure sixe
set(gcf, 'Position', [0 0 2500 700])


