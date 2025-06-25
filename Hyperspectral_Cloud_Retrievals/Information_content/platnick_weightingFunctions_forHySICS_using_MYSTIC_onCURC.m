%% Compute weighting functions using libRadtran's MYSTIC Monte Carlo Model

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
        'Hyperspectral_Cloud_Retrievals/HySICS/Monte_Carlo/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/testing_MYSTIC/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Monte_Carlo/'];


    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/testing_MYSTIC/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    inputs.folderpath_reflectance = '/scratch/alpine/anbu8374/Weighting_functions/';



    % Define the folder path where all .INP files will be saved
    inputs.folderpath_inp = '/scratch/alpine/anbu8374/Weighting_functions/';

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




%% set up the inputs to create an INP file for MYSTIC!


[inputs, spec_response] = create_uvSpec_MYSTIC_inputs_for_HySICS(inputs);


%% Write the INP files



tic

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % length of each independent variable
    num_rEff = length(inputs.RT.re);
    num_tauC = length(inputs.RT.tau_c);
    num_wl = size(inputs.RT.wavelengths2run,1);

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
            inputs.RT.compute_weighting_functions, inputs.which_computer, idx_unique_indVars(nn));

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
    % num_rTop = length(inputs.RT.r_top);
    % num_rBot = length(inputs.RT.r_bot);
    % num_tauC = length(inputs.RT.tau_c);
    num_wl = size(inputs.RT.wavelengths2run,1);
    num_tau_layers = inputs.RT.n_layers;

    % num_INP_files = num_rTop*num_rBot*num_tauC*num_wl;
    num_INP_files = num_wl*num_tau_layers;


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

    % the wc_filenames should be the same for different wavelengths
    wc_filename = repmat(wc_filename, num_wl, 1);
    

    % Compute the optical depth of each layer
    inputs.RT.tau_layers = (lwc.*ext_bulk_coeff_per_LWC.*(inputs.RT.z_edges(2) - inputs.RT.z_edges(1)))';  % the optical depth of each layer, starting from cloud top



    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for rt = 1:num_rTop
    %   for rb = 1:num_rBot
    %       for tt = 1:num_tauC
    %           for ww = 1:num_wl
    % changing_variables = [reshape(repmat(inputs.RT.r_top, num_rBot*num_tauC*num_wl,1), [],1),...
    %     repmat(reshape(repmat(inputs.RT.r_bot, num_tauC*num_wl,1), [],1), num_rTop, 1),...
    %     repmat(reshape(repmat(inputs.RT.tau_c, num_wl,1), [],1), num_rBot*num_rTop, 1),...
    %     repmat(inputs.RT.wavelengths2run, num_rTop*num_rBot*num_tauC, 1)];

    changing_variables = [reshape(repmat(inputs.RT.wavelengths2run(:,1)', num_tau_layers,1), [],1),...
        repmat(flipud(cumsum(fliplr(inputs.RT.tau_layers))'), num_wl,1)];


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, reshape(repmat((1:num_wl), num_tau_layers, 1), [],1)];

    
    % define the basename
    mc_basename = cell(num_INP_files, 1);



    % Now write all the INP files
    parfor nn = 1:num_INP_files
        %     for nn = 1:num_INP_files

        % update the monte carlo basename
        mc_basename{nn} = [inputs.folderpath_inp, 'mc_',num2str(inputs.RT.wavelengths2run(1)),...
            'nm_layers1-', num2str(num_INP_files - (nn-1))];

        % set the wavelengths for each file
        wavelengths = changing_variables(nn);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = ['monteCarlo_',num2str(mean(wavelengths)), '_','nm_rTop_', num2str(inputs.RT.r_top),...
            '_rBot_', num2str(inputs.RT.r_bot),'_tauC_', num2str(round(changing_variables(nn,2),4)),...
            '_layers1-', num2str(num_INP_files - (nn-1)), '.INP'];


        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];




        % ------------------ Write the INP File --------------------
        write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn}, mc_basename{nn});


    end



end

toc



%% Run the MYSTIC program

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
    [ds,~,~] = readUVSPEC_ver2(inputs.folderpath_inp, mc_basename{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);


    % Store the reflectance
    Refl_model(nn, :) = ds.reflectivity;       % reflectance is in units of 1/sr

    % *** Monte Carlo Cannot compute azimuthally avg. reflectance ***
    % Store the azimuthally averaged reflectance
    % Azimuthal average at this particular zenith angle
    % Refl_model(nn, :) = ds.reflectivity.az_avg;       % reflectance is in units of 1/sr

    % % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    % [Refl_model(nn), ~] = reflectanceFunction(inputSettings(2,:), ds,...
    %     spec_response_value(changing_variables(nn,end),:));



end






toc




%% Compute the weighting functions using Platnick (2000)
% Equation 4 defines the weighting function as the normalized derivative of
% reflectance with respect to optical depth.


% reshape Refl_model
Refl_model = reshape(Refl_model, inputs.RT.n_layers, size(inputs.RT.wavelengths2run,1));




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
        'Hyperspectral_Cloud_Retrievals/HySICS/Monte_Carlo/'];



elseif strcmp(inputs.which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Monte_Carlo/'];


elseif strcmp(inputs.which_computer,'curc')==true

    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

    inputs.folderpath_2save = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/Monte_Carlo/';



end


% If the folder path doesn't exit, create a new directory
if ~exist(inputs.folderpath_2save, 'dir')

    mkdir(inputs.folderpath_2save)

end



rev = 1;

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous')==true


    filename = [inputs.folderpath_2save,'monteCarlo_HySICS_reflectance_for_weightingFunctions_',...
            num2str(round(inputs.RT.wavelengths2run(1))), 'nm_inhomogeneous_droplet_profile_sim-ran-on-',...
            char(datetime("today")), '_rev', num2str(rev),'.mat'];


else

    filename = [inputs.folderpath_2save,'monteCarlo_HySICS_reflectance_for_weightingFunctions_',...
            num2str(round(inputs.RT.wavelengths2run(1))), 'nm_homogeneous_droplet_profile_sim-ran-on-',...
            char(datetime("today")), '_rev', num2str(rev),'.mat'];

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

