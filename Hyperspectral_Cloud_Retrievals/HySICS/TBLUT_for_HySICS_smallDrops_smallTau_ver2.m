%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using simulated HySICS data


% *** CURRENT FORWARD MODEL UNCERTAINTIES CONSIDERED ***
% (1) Adiabatic droplet profile assumption
% (2) Cloud top height assumption


% ---------------- INPUTS ---------------
% ---------------------------------------
% (1) simulated_reflectance - HySICS simulated measurements

% (2) folderpaths - this is a structure with all the folder paths needed to
% read and store files, including the path to the simulated data set, the
% location to save the retirevals, and the location to save the INP files

% (3) print_status_updates - true or false that tells the function to
% display updates on the progress of the calculation

% (4) print_libRadtran_err - true or false that tells the function to
% write and save the libRadtran error message file




% By Andrew John Buggee

%%

function tblut_retrieval_2 = TBLUT_for_HySICS_smallDrops_smallTau_ver2(simulated_measurements, folder_paths,...
    tblut_retrieval_1, print_status_updates, print_libRadtran_err)


%% unpack folder_paths

which_computer = folder_paths.which_computer;


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_tblut_2 = create_HySICS_inputs_TBLUT_smallTau_smallZenithAngles(simulated_measurements.inputs, tblut_retrieval_1,...
    print_libRadtran_err);


%% Update the cloud top height!
% ** VOCALS-REx in-situ measurements result in a mean cloud top height
% of 1203 meters and a mean cloud depth of about 230 meters
% ** testing the retrieval when I lack knowledge of cloud top precisely **

% load the set of VOCALS-REx in-situ observations
if strcmp(which_computer, 'anbu8374')==true

    cloud_top_obs = load(['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

elseif strcmp(which_computer, 'andrewbuggee')==true

    cloud_top_obs = load(['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Presentations_and_Papers/paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

elseif strcmp(which_computer, 'curc')==true

    cloud_top_obs = load(['/projects/anbu8374/Matlab-Research/Presentations_and_Papers/',...
        'paper_2/VR_cloud_top_height_obs_19-Jan-2026.mat']);

end

cth_mean = mean(cloud_top_obs.cloudTopHeight)/1e3;                  % km
cld_depth_mean = mean(cloud_top_obs.cloudDepth)/1e3;                % km

inputs_tblut_2.RT.z_topBottom = [cth_mean, (cth_mean - cld_depth_mean)];         % kilometers

% update depenent variables
inputs_tblut_2.RT.cloud_depth = cld_depth_mean;                % kilometers
inputs_tblut_2.RT.H = inputs_tblut_2.RT.z_topBottom(1) - inputs_tblut_2.RT.z_topBottom(2);                                % km - geometric thickness of cloud



%% Find the measurements closest to the bands to run

[~, idx_1] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut_2.bands2run(1)));

[~, idx_2] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut_2.bands2run(2)));

[~, idx_3] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut_2.bands2run(3)));

[~, idx_4] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut_2.bands2run(4)));

% error if the values found are at least 15nm from the intended wavlengths
if abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_1,:)) - 500)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 500 nm', newline])

elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_2,:)) - 1680)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 1680 nm', newline])

elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_3,:)) - 2131)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 2131 nm', newline])

elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_4,:)) - 2226)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 2226 nm', newline])

else

    % Then we set the bands to run to be to ones found to be closest to the
    % desired bands out of the measurement bands provided
    inputs_tblut_2.bands2run_from_set_of_measurements = [idx_1, idx_2, idx_3, idx_4];
    inputs_tblut_2.bands2plot = inputs_tblut_2.bands2run;

    % ---- Define the wavelengths ----
    inputs_tblut_2.RT.wavelengths2run = simulated_measurements.inputs.RT.wavelengths2run(inputs_tblut_2.bands2run_from_set_of_measurements,:);



end
%% ----- Create .INP files for HySICS TBLUT -----


if inputs_tblut_2.flags.writeINPfiles == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % *** This time we dont' sweep through the values, we compute the
    % reflectance at each as if it were a state vector

    % length of each independent variable
    num_rEff = length(inputs_tblut_2.RT.re);
    % num_tauC = length(inputs_tblut_2.RT.tau_c);
    num_wl = length(inputs_tblut_2.bands2run);

    num_INP_files = num_rEff*num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for xx = 1:num_states
    %       for ww = 1:num_wl
    changing_variables = [reshape(repmat(inputs_tblut_2.RT.re, num_wl, 1), [], 1),...
        reshape(repmat(inputs_tblut_2.RT.tau_c, num_wl, 1), [], 1),...
        repmat(inputs_tblut_2.RT.wavelengths2run, num_rEff, 1)];


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_rEff, 1)];


    % First, write all the wc files
    temp_names = cell(num_rEff, 1);
    wc_filename = cell(num_INP_files, 1);

    % Define the water cloud folder path location
    wc_folder_path = folder_paths.libRadtran_water_cloud_files;

    % Define the mie folder path
    mie_folder_path = folder_paths.libRadtran_mie_folder;

    % only jump on indexes where there is a unique r and tau pair
    idx_parFor = 1:(num_rEff-1):num_INP_files;

    parfor nn = 1:numel(idx_parFor)
        % for nn = 1:num_rEff

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------


        temp = write_wc_file(changing_variables(idx_parFor(nn), 1), changing_variables(idx_parFor(nn),2),...
            inputs_tblut_2.RT.z_topBottom,inputs_tblut_2.RT.lambda_forTau, inputs_tblut_2.RT.distribution_str,...
            inputs_tblut_2.RT.distribution_var,inputs_tblut_2.RT.vert_homogeneous_str, inputs_tblut_2.RT.parameterization_str,...
            inputs_tblut_2.RT.indVar, inputs_tblut_2.compute_weighting_functions, which_computer, nn+(nn-1), 1,...
            wc_folder_path, mie_folder_path);

        temp_names{nn} = temp{1};

    end

    % set the odd and even values to have the same file names
    for nn = 1:num_rEff
        wc_filename(idx_parFor(nn):idx_parFor(nn)+(num_rEff-1)) = temp_names(nn);
    end



    % define the INP folder location
    inp_folder_path = folder_paths.libRadtran_inp;

    % define the libRadtran data path
    libRadtran_data_path = folder_paths.libRadtran_data;



    % Now write all the INP files
    parfor nn = 1:num_INP_files
        % for nn = 1:num_INP_files


        % set the wavelengths for each file
        wavelengths = changing_variables(nn, end-2:end-1);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = ['TBLUT_retrieval_', num2str(mean(wavelengths)), '_',...
            'nm_rEff_', num2str(changing_variables(nn,1)), '_tauC_', num2str(changing_variables(nn,2)), '_',...
            inputs_tblut_2.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(inp_folder_path, libRadtran_data_path, wc_folder_path, inputFileName{nn},...
            inputs_tblut_2, wavelengths, wc_filename{nn});


    end




else

    % if the files already exist, just grab the names!
    error([newline, 'Dont know how to do this!', newline])


end







%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs_tblut_2.flags.runUVSPEC == true


    % Read the solar flux file over the wavelength range specified
    wavelength_vec = [min(inputs_tblut_2.RT.wavelengths2run,[],"all"), max(inputs_tblut_2.RT.wavelengths2run, [], "all")];

    [source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs_tblut_2.RT.source_file);   % W/nm/m^2

    % we will add and subtract a small fraction of the source file resolution
    % to ensure rounding errors don't cause an issue when selecting the
    % wavelengths needed from the source file
    wl_perturb = inputs_tblut_2.RT.source_file_resolution/3;   % nm


    % Let's only keep the source flux values needed for the parfor loop
    % for the first channel
    idx_wl_1 = source_wavelength>=(changing_variables(1,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables(1,end-1) + wl_perturb);

    % for the second channel
    idx_wl_2 = source_wavelength>=(changing_variables(2,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables(2,end-1) + wl_perturb);

    % for the third channel
    idx_wl_3 = source_wavelength>=(changing_variables(3,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables(3,end-1) + wl_perturb);

    % for the fourth channel
    idx_wl_4 = source_wavelength>=(changing_variables(4,end-2) - wl_perturb) &...
        source_wavelength<=(changing_variables(4,end-1) + wl_perturb);

    source_flux_for_parForloop = [source_flux(idx_wl_1), source_flux(idx_wl_2),...
        source_flux(idx_wl_3), source_flux(idx_wl_4)]; % the only source flux values needed for the loop



    spec_response = simulated_measurements.spec_response.value;

    % Only keep the spectral response functions needed in the parfor loop
    spec_response = spec_response(inputs_tblut_2.bands2run_from_set_of_measurements, :);


    % store the reflectances
    Refl_model_tblut_2 = zeros(num_INP_files, 1);


    if print_status_updates==true

        parfor nn = 1:num_INP_files
            % for nn = 1:num_INP_files


            disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(inp_folder_path, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(inp_folder_path, outputFileName{nn}, inputs_tblut_2,...
                inputs_tblut_2.RT.compute_reflectivity_uvSpec);


            % compute the reflectance
            [Refl_model_tblut_2(nn), ~] = reflectanceFunction_ver2(inputs_tblut_2, ds,...
                source_flux_for_parForloop(:,changing_variables(nn,end)), spec_response(changing_variables(nn,end),:)');



        end


    else


        parfor nn = 1:num_INP_files
            % for ww = 1:num_INP_files


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(inp_folder_path, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(inp_folder_path, outputFileName{nn}, inputs_tblut_2,...
                inputs_tblut_2.RT.compute_reflectivity_uvSpec);


            % compute the reflectance
            [Refl_model_tblut_2(nn), ~] = reflectanceFunction_ver2(inputs_tblut_2, ds,...
                source_flux_for_parForloop(:,changing_variables(nn,end)), spec_response(changing_variables(nn,end),:)');



        end


    end










    % save the calculated reflectances and the inputs
    save(folder_paths.saveOutput_filename, "inputs_tblut_2", "Refl_model_tblut_2", '-append'); % save inputSettings to the same folder as the input and output file



elseif inputs_tblut_2.flags.runUVSPEC == false

    load([inputs_tblut_2.savedCalculations_folderName,inputs_tblut_2.saveCalculations_fileName] ,'inputs_tblut_2','Refl_model_tblut_2');

end


%% Let's comapre the newly calculated reflectances with the measured reflectances

% Grab observations
if isfield(simulated_measurements, 'Refl_model_with_noise')==true

    % for band 1...
    observations_band1 = simulated_measurements.Refl_model_with_noise(inputs_tblut_2.bands2run_from_set_of_measurements(1));
    % for band 2...
    observations_band2 = simulated_measurements.Refl_model_with_noise(inputs_tblut_2.bands2run_from_set_of_measurements(2));
    % for band 3...
    observations_band3 = simulated_measurements.Refl_model_with_noise(inputs_tblut_2.bands2run_from_set_of_measurements(3));
    % for band 4...
    observations_band4 = simulated_measurements.Refl_model_with_noise(inputs_tblut_2.bands2run_from_set_of_measurements(4));

elseif isfield(simulated_measurements, 'Refl_model')==true

    % for band 1...
    observations_band1 = simulated_measurements.Refl_model(inputs_tblut_2.bands2run_from_set_of_measurements(1));
    % for band 2...
    observations_band2 = simulated_measurements.Refl_model(inputs_tblut_2.bands2run_from_set_of_measurements(2));
    % for band 3...
    observations_band3 = simulated_measurements.Refl_model(inputs_tblut_2.bands2run_from_set_of_measurements(3));
    % for band 4...
    observations_band4 = simulated_measurements.Refl_model(inputs_tblut_2.bands2run_from_set_of_measurements(4));

elseif isfield(simulated_measurements, 'Refl_model_2save')==true

    % for band 1...
    observations_band1 = simulated_measurements.Refl_model_2save(inputs_tblut_2.bands2run_from_set_of_measurements(1));
    % for band 2...
    observations_band2 = simulated_measurements.Refl_model_2save(inputs_tblut_2.bands2run_from_set_of_measurements(2));
    % for band 3...
    observations_band3 = simulated_measurements.Refl_model_2save(inputs_tblut_2.bands2run_from_set_of_measurements(3));
    % for band 4...
    observations_band4 = simulated_measurements.Refl_model_2save(inputs_tblut_2.bands2run_from_set_of_measurements(4));


end

% take the difference between what was calculated and what was measured
% reshape the libRadTran forward modeled reflectance to be an array where one
% dimension varies with optical depth and another varies with
% effective radius
modelRefl_band1_array = Refl_model_tblut_2(1:num_wl:length(Refl_model_tblut_2));
modelRefl_band2_array = Refl_model_tblut_2(2:num_wl:length(Refl_model_tblut_2));
modelRefl_band3_array = Refl_model_tblut_2(3:num_wl:length(Refl_model_tblut_2));
modelRefl_band4_array = Refl_model_tblut_2(4:num_wl:length(Refl_model_tblut_2));

% now take the difference between the true measurements and the
% calculations and find the minimum RMS value
tblut_retrieval_2.minRMS = (mean( [(modelRefl_band1_array - observations_band1).^2 ,...
    (modelRefl_band2_array - observations_band2).^2,...
    (modelRefl_band3_array - observations_band3).^2,...
    (modelRefl_band4_array - observations_band4).^2], 2 )).^(0.5);

tblut_retrieval_2.states_2check = [inputs_tblut_2.RT.re', inputs_tblut_2.RT.tau_c'];







end
