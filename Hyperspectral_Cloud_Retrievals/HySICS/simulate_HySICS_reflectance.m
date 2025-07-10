%% Create simulated HySICS measurements
% Include Gaussian random noise


% By Andrew John Buggee


clear variables


%% Define the path location where INP files will be stored, and where Reflectances will be stored

clear inputs

% Determine which computer this is being run on
inputs.which_computer = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(inputs.which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];


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
    inputs.folderpath_reflectance = '/scratch/alpine/anbu8374/hyperspectral_retrieval/HySICS/';



    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = '/scratch/alpine/anbu8374/hyperspectral_retrieval/HySICS/INP-OUT/';

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



%% Do you want to model 2 parameters for the droplet profile (r_top and r_bot)
% or 3 (r_top, r_middle, r_bot)

inputs.RT.num_re_parameters = 2;

%% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file

[inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, false, [], 'exact');

inputs.calc_type = 'simulated_spectra';

%% Set the total column water vapor?

inputs.RT.modify_total_columnWaterVapor = true;             % modify the full column
inputs.RT.waterVapor_column = 20;    % mm

inputs.RT.modify_aboveCloud_columnWaterVapor = false;         % don't modify the column above the cloud

% determine the above cloud column water vapor amount
inputs.RT.aboveCloud_cwv = aboveCloud_CWV_simulated_hysics_spectra(inputs);


%% Write each INP file



% num wavelengths
num_wl = length(inputs.bands2run);



idx = 0;


tic

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % length of each independent variable
    num_rEff = length(inputs.RT.re);
    num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);

    num_INP_files = num_rEff*num_tauC*num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for rr = 1:num_rEff
    %   for tt = 1:num_tauC
    %       for ww = 1:num_wl
    changing_variables = [reshape(repmat(inputs.RT.re, num_tauC*num_wl,1), [],1),...
        repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rEff, 1),...
        repmat(inputs.RT.wavelengths2run, num_rEff*num_tauC, 1)];

    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_rEff*num_tauC, 1)];


    % First, write all the wc files
    temp_names = cell(num_rEff*num_tauC, 1);
    wc_filename = cell(num_INP_files, 1);

    % only jump on indexes where there is a unique r and tau pair
    idx_unique_indVars = 1:num_wl:num_INP_files;

    parfor nn = 1:length(idx_unique_indVars)

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------

        temp = write_wc_file(changing_variables(idx_unique_indVars(nn), 1),...
            changing_variables(idx_unique_indVars(nn), 2),inputs.RT.z_topBottom,...
            inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
            inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
            inputs.RT.indVar, inputs.compute_weighting_functions, inputs.which_computer,...
            idx_unique_indVars(nn), inputs.RT.num_re_parameters);

        temp_names{nn} = temp{1};

    end

    % the wc_filenames should be the same for different wavelengths
    for ww = 0:num_wl-1
        wc_filename(idx_unique_indVars+ww) = temp_names;
    end





    % Now write all the INP files
    parfor nn = 1:num_INP_files


        % set the wavelengths for each file
        wavelengths = changing_variables(nn, end-2:end-1);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = [num2str(mean(wavelengths)), '_',...
            'nm_rEff_', num2str(changing_variables(nn,1)), '_tauC_', num2str(changing_variables(nn,2)), '_',...
            inputs.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn});


    end







elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------



    % length of each independent variable
    % length of each independent variable
    num_rTop = length(inputs.RT.r_top);
    num_rBot = length(inputs.RT.r_bot);
    num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);

    num_INP_files = num_rTop*num_rBot*num_tauC*num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    if inputs.RT.num_re_parameters==2
        % changing variable steps through reff, tauC, and wavelength
        % in for loop speak, it would be:
        % for rt = 1:num_rTop
        %   for rb = 1:num_rBot
        %       for tt = 1:num_tauC
        %           for ww = 1:num_wl
        changing_variables = [reshape(repmat(inputs.RT.r_top, num_rBot*num_tauC*num_wl,1), [],1),...
            repmat(reshape(repmat(inputs.RT.r_bot, num_tauC*num_wl,1), [],1), num_rTop, 1),...
            repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rBot*num_rTop, 1),...
            repmat(inputs.RT.wavelengths2run, num_rTop*num_rBot*num_tauC, 1)];

    elseif inputs.RT.num_re_parameters==3

        % changing variable steps through reff, tauC, and wavelength
        % in for loop speak, it would be:
        % for rt = 1:num_rTop
        %   for rm = 1:num_rMiddle
        %       for rb = 1:num_rBot
        %           for tt = 1:num_tauC
        %               for ww = 1:num_wl
        changing_variables = [reshape(repmat(inputs.RT.r_top, num_rMiddle*num_rBot*num_tauC*num_wl,1), [],1),...
            repmat(reshape(repmat(inputs.RT.r_bot, num_tauC*num_wl,1), [],1), num_rTop, 1),...
            repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rBot*num_rTop, 1),...
            repmat(inputs.RT.wavelengths2run, num_rTop*num_rBot*num_tauC, 1)];

    else

        error([newline, 'How many free parameters are used to model the droplet profile?', newline])

    end


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_rTop*num_rBot*num_tauC, 1)];

    % First, write all the wc files
    temp_names = cell(num_rTop*num_rBot*num_tauC, 1);
    wc_filename = cell(num_INP_files, 1);

    % only jump on indexes where there is a unique r and tau pair
    idx_unique_indVars = 1:num_wl:num_INP_files;

    parfor nn = 1:length(idx_unique_indVars)

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------

        re = create_droplet_profile2([changing_variables(idx_unique_indVars(nn), 1),...
            changing_variables(idx_unique_indVars(nn), 2)],...
            inputs.RT.z, inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector


        temp = write_wc_file(re, changing_variables(idx_unique_indVars(nn), 3),...
            inputs.RT.z_topBottom,inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
            inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
            inputs.RT.indVar, inputs.compute_weighting_functions, inputs.which_computer,...
            idx_unique_indVars(nn), inputs.RT.num_re_parameters);

        temp_names{nn} = temp{1};

    end

    % the wc_filenames should be the same for different wavelengths
    for ww = 0:num_wl-1
        wc_filename(idx_unique_indVars+ww) = temp_names;
    end


    % Now write all the INP files
    if inputs.RT.num_re_parameters==2

        parfor nn = 1:num_INP_files
            % for nn = 1:num_INP_files


            % set the wavelengths for each file
            wavelengths = changing_variables(nn, end-2:end-1);

            % ------------------------------------------------
            % ---- Define the input and output filenames! ----
            % ------------------------------------------------
            % input_names need a unique identifier. Let's give them the nn value so
            % they can be traced, and are writen over in memory


            inputFileName{nn} = [num2str(mean(wavelengths)), '_','nm_rTop_', num2str(changing_variables(nn,1)),...
                '_rBot_', num2str(changing_variables(nn,2)),'_tauC_', num2str(changing_variables(nn,3)), '_',...
                inputs.RT.atm_file(1:end-4),'.INP'];



            outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


            % ------------------ Write the INP File --------------------
            write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
                wavelengths, wc_filename{nn});


        end




    elseif inputs.RT.num_re_parameters==3



        parfor nn = 1:num_INP_files
            % for nn = 1:num_INP_files


            % set the wavelengths for each file
            wavelengths = changing_variables(nn, end-2:end-1);

            % ------------------------------------------------
            % ---- Define the input and output filenames! ----
            % ------------------------------------------------
            % input_names need a unique identifier. Let's give them the nn value so
            % they can be traced, and are writen over in memory


            inputFileName{nn} = [num2str(mean(wavelengths)), '_','nm_rTop_', num2str(changing_variables(nn,1)),...
                '_rMiddle_', num2str(changing_variables(nn,2)), '_rBot_', num2str(changing_variables(nn,3)),...
                '_tauC_', num2str(changing_variables(nn,4)), '_', inputs.RT.atm_file(1:end-4),'.INP'];



            outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


            % ------------------ Write the INP File --------------------
            write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
                wavelengths, wc_filename{nn});


        end

    else

        error([newline, 'How many free parameters are used to model the droplet profile?', newline])

    end



end

toc






%% Calculate Reflectance

% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(inputs.RT.wavelengths2run,[],"all"), max(inputs.RT.wavelengths2run, [], "all")];
[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = inputs.RT.source_file_resolution/3;   % nm



% define only the spec_response so the wavelengths are passed into the
% memory of the parallel for loop
spec_response_value = spec_response.value;


tic



% store the reflectances
Refl_model = zeros(num_INP_files, 1);


parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


    disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------

    % compute INP file
    runUVSPEC_ver2(inputs.folderpath_inp, inputFileName{nn}, outputFileName{nn},...
        inputs.which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(inputs.folderpath_inp, outputFileName{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(changing_variables(nn,4) - wl_perturb) &...
        source_wavelength<=(changing_variables(nn,5) + wl_perturb);


    % Store the Radiance
    %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    [Refl_model(nn), ~] = reflectanceFunction_ver2(inputs, ds,...
        source_flux(idx_wl), spec_response_value(changing_variables(nn,end),:));



end






toc


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
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];



elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];


elseif strcmp(inputs.which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    inputs.folderpath_2save = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/';



end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_2save, 'dir')

    mkdir(inputs.folderpath_2save)

end



rev = 1;

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true

    if strcmp(inputs.calc_type, 'forward_model_calcs_forRetrieval')==true

        filename = [inputs.folderpath_2save,'forward_model_calcs_forRetrieval_HySICS_reflectance_',...
            'inhomogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

    elseif strcmp(inputs.calc_type, 'simulated_spectra')==true

        filename = [inputs.folderpath_2save,'simulated_measurement_HySICS_reflectance_',...
            'inhomogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

    end


else

    if strcmp(inputs.calc_type, 'forward_model_calcs_forRetrieval')==true

        filename = [inputs.folderpath_2save,'forward_model_calcs_forRetrieval_HySICS_reflectance_',...
            'homogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

    elseif strcmp(inputs.calc_type, 'simulated_spectra')==true

        filename = [inputs.folderpath_2save,'simulated_measurement_HySICS_reflectance_',...
            'homogeneous_droplet_profile_sim-ran-on-',char(datetime("today")), '_rev', num2str(rev),'.mat'];

    end

end


while isfile(filename)
    rev = rev+1;
    if rev<10
        filename = [filename(1:end-5), num2str(rev),'.mat'];
    elseif rev>10
        filename = [filename(1:end-6), num2str(rev),'.mat'];
    end
end


save(filename, "Refl_model","inputs", "spec_response", "changing_variables");



%% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
% inputs.measurement.uncert = [0.003, 0.01:0.01:0.1];
inputs.measurement.uncert = 0;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements.
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between
% +/- 5% of the measured value

% To sample a normal distribution with mean mu, and standard deviation s,
% we compute the following: y = s * randn() + mu

% We define the standard deviation as the measurement uncertainty divided
% by three. Therefore, after sample a large number of times, 99% of
% measurements will be within +/- measurement uncertainy of the mean

if any(inputs.measurement.uncert > 0)

    inputs.measurement.standard_dev = inputs.measurement.uncert/3;       % this is still just a fraction

    for uu = 1:length(inputs.measurement.uncert)

        clear Refl_model_with_noise Refl_model_uncert


        Refl_model_with_noise = (inputs.measurement.standard_dev(uu) .* Refl_model) .* randn(size(Refl_model))...
            + Refl_model;


        % define the synthetic relfectance uncertainty
        Refl_model_uncert = inputs.measurement.uncert(uu) .* Refl_model_with_noise;    % 1/sr



        % Save the output with added measurement uncertainty
        ver = 1;

        if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true

            if strcmp(inputs.calc_type, 'forward_model_calcs_forRetrieval')==true

                filename = [inputs.folderpath_2save,'forward_model_calcs_forRetrieval_HySICS_reflectance_',...
                    'inhomogeneous_droplet_profile_with_',num2str(inputs.measurement.uncert(uu)*100),...
                    '%_uncertainty_sim-ran-on-', char(datetime("today")), '_rev', num2str(ver),'.mat'];

            elseif strcmp(inputs.calc_type, 'simulated_measurement')==true

                filename = [inputs.folderpath_2save,'simulated_HySICS_reflectance_inhomogeneous_droplet_profile_',...
                    'with_',num2str(inputs.measurement.uncert(uu)*100), '%_uncertainty_sim-ran-on-',...
                    char(datetime("today")), '_rev', num2str(ver),'.mat'];

            end



        else

            if strcmp(inputs.calc_type, 'forward_model_calcs_forRetrieval')==true

                filename = [inputs.folderpath_2save,'forward_model_calcs_forRetrieval_HySICS_reflectance_',...
                    'homogeneous_droplet_profile_with_',num2str(inputs.measurement.uncert(uu)*100),...
                    '%_uncertainty_sim-ran-on-', char(datetime("today")), '_rev', num2str(ver),'.mat'];

            elseif strcmp(inputs.calc_type, 'simulated_measurement')==true

                filename = [inputs.folderpath_2save,'simulated_HySICS_reflectance_homogeneous_droplet_profile_',...
                    'with_',num2str(inputs.measurement.uncert(uu)*100), '%_uncertainty_sim-ran-on-',...
                    char(datetime("today")), '_rev', num2str(ver),'.mat'];

            end


        end


        while isfile(filename)
            ver = ver+1;
            filename = [filename(1:end-5), num2str(ver), '.mat'];
        end


        save(filename, "Refl_model_with_noise", "Refl_model_uncert","inputs",...
            "changing_variables", "spec_response");

    end

end






%% Plot the results


figure;
if size(inputs.RT.wavelengths2run,1)>1 && size(inputs.RT.wavelengths2run,2)>1

    if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true && ...
            num_rTop==1 && num_rBot==1 && num_tauC==1

        % There is one state vector computed for a range of wavelengths
        plot(mean(inputs.RT.wavelengths2run,2),Refl_model,...
            '.-', 'linewidth', 1, 'markersize', 35, 'Color', mySavedColors(1, 'fixed'))


        title('Simulated Reflectance - liquid water cloud','Interpreter', 'latex')
        subtitle(['$r_{top} = $',num2str(round(inputs.RT.r_top,1)), ' $\mu m$, $r_{bot} = $',...
            num2str(round(inputs.RT.r_bot,1)), ' $\mu m$, $\tau_c = $', num2str(round(inputs.RT.tau_c,1))],...
            'Interpreter', 'latex')



    elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true && ...
            num_rTop>1 && num_rBot>1 && num_tauC>1

        % indexes to plot
        r_top_2Plot = changing_variables(:,1)==9;
        r_bot_2Plot = changing_variables(:,2)==4;

        for nn = 1:size(Refl_model,4)

            plot(mean(inputs.RT.wavelengths2run,2),...
                reshape(Refl_model(:,r_top_2Plot, r_bot_2Plot, nn), 1, []),...
                '.-', 'linewidth', 1, 'markersize', 35, 'Color', mySavedColors(nn, 'fixed'))

            hold on

            legend_str{nn} = ['$\tau_c = $', num2str(inputs.RT.tau_c(nn))];

        end



    elseif size(Refl_model,2)==1 && size(Refl_model,3)==1 && size(Refl_model,4)>1


        for nn = 1:size(Refl_model,4)

            plot(mean(inputs.RT.wavelengths2run,2),...
                reshape(Refl_model(:,1, 1, nn), length(inputs.bands2run), []),...
                '.-', 'linewidth', 1, 'markersize', 35, 'Color', mySavedColors(nn, 'fixed'))

            hold on

            legend_str{nn} = ['$\tau_c = $', num2str(inputs.RT.tau_c(nn))];

        end



    end

end

hold on

grid on; grid minor
xlabel('Wavelength (nm)','Interpreter', 'latex')
ylabel('Reflectance (1/sr)','Interpreter', 'latex')
set(gcf, 'Position', [0 0 1000 1000])
%legend(legend_str,  'Interpreter', 'latex', 'Fontsize', 30', 'location', 'best')
% title(['Simulated Reflectance - liquid water cloud - $r_e = $', num2str(re_2plot), ' $\mu m$, $\tau_c = $',...
%     num2str(tau_2plot)], 'Interpreter', 'latex')
