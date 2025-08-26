%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using simulated HySICS data


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

function tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_measurements, folder_paths, print_status_updates, print_libRadtran_err)


%% unpack folder_paths

which_computer = folder_paths.which_computer;


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_tblut = create_HySICS_inputs_TBLUT(folder_paths, simulated_measurements.inputs, print_libRadtran_err);



%% Find the measurements closest to the bands to run

[~, idx_1] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut.bands2run(1)));

[~, idx_2] = min(abs(simulated_measurements.inputs.bands2run - inputs_tblut.bands2run(2)));

% error if the values found are at least 15nm from the intended wavlengths
if abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_1,:)) - 502)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 650 nm', newline])

elseif abs(mean(simulated_measurements.inputs.RT.wavelengths2run(idx_2,:)) - 2130)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 650 nm', newline])

else

    % Then we set the bands to run to be to ones found to be closest to the
    % desired bands out of the measurement bands provided
    inputs_tblut.bands2run_from_set_of_measurements = [idx_1, idx_2];
    inputs_tblut.bands2plot = inputs_tblut.bands2run;

    % ---- Define the wavelengths ----
    inputs_tblut.RT.wavelengths2run = simulated_measurements.inputs.RT.wavelengths2run(inputs_tblut.bands2run_from_set_of_measurements,:);



end
%% ----- Create .INP files for HySICS TBLUT -----


if inputs_tblut.flags.writeINPfiles == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % length of each independent variable
    num_rEff = length(inputs_tblut.RT.re);
    num_tauC = length(inputs_tblut.RT.tau_c);
    num_wl = length(inputs_tblut.bands2run);

    num_INP_files = num_rEff*num_tauC*num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through reff, tauC, and wavelength
    % in for loop speak, it would be:
    % for rr = 1:num_rEff
    %   for tt = 1:num_tauC
    %       for ww = 1:num_wl
    changing_variables = [reshape(repmat(inputs_tblut.RT.re, num_tauC*num_wl,1), [],1),...
        repmat(reshape(repmat(inputs_tblut.RT.tau_c, num_wl,1), [],1), num_rEff, 1),...
        repmat(inputs_tblut.RT.wavelengths2run, num_rEff*num_tauC, 1)];

    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:num_wl)', num_rEff*num_tauC, 1)];


    % First, write all the wc files
    temp_names = cell(num_rEff*num_tauC, 1);
    wc_filename = cell(num_INP_files, 1);

    % Define the water cloud folder path location
    wc_folder_path = folder_paths.libRadtran_water_cloud_files;

    % Define the mie folder path
    mie_folder_path = folder_paths.libRadtran_mie_folder;

    % only jump on indexes where there is a unique r and tau pair

    parfor nn = 1:num_rEff*num_tauC
        % for nn = 1:num_rEff*num_tauC

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------


        temp = write_wc_file(changing_variables(2*nn -1, 1), changing_variables(2*nn -1,2),...
            inputs_tblut.RT.z_topBottom,inputs_tblut.RT.lambda_forTau, inputs_tblut.RT.distribution_str,...
            inputs_tblut.RT.distribution_var,inputs_tblut.RT.vert_homogeneous_str, inputs_tblut.RT.parameterization_str,...
            inputs_tblut.RT.indVar, inputs_tblut.compute_weighting_functions, which_computer, nn+(nn-1), 1,...
            wc_folder_path, mie_folder_path);

        temp_names{nn} = temp{1};

    end

    % set the odd and even values to have the same file names
    wc_filename(1:num_wl:num_INP_files, 1) = temp_names;
    wc_filename(2:num_wl:num_INP_files, 1) = temp_names;


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
            inputs_tblut.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(inp_folder_path, libRadtran_data_path, wc_folder_path, inputFileName{nn},...
            inputs_tblut, wavelengths, wc_filename{nn});


    end




else

    % if the files already exist, just grab the names!
    error([newline, 'Dont know how to do this!', newline])


end







%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs_tblut.flags.runUVSPEC == true


    % Read the solar flux file over the wavelength range specified
    wavelength_vec = [min(inputs_tblut.RT.wavelengths2run,[],"all"), max(inputs_tblut.RT.wavelengths2run, [], "all")];

    [source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs_tblut.RT.source_file);   % W/nm/m^2

    % we will add and subtract a small fraction of the source file resolution
    % to ensure rounding errors don't cause an issue when selecting the
    % wavelengths needed from the source file
    wl_perturb = inputs_tblut.RT.source_file_resolution/3;   % nm


    % Let's only keep the source flux values needed for the parfor loop




    spec_response = simulated_measurements.spec_response.value;

    % Only keep the spectral response functions needed in the parfor loop
    spec_response = spec_response(inputs_tblut.bands2run_from_set_of_measurements, :);


    % store the reflectances
    Refl_model_tblut = zeros(num_INP_files, 1);


    if print_status_updates==true

        parfor nn = 1:num_INP_files
            % for ww = 1:size(inputs.RT.wavelengths2run, 1)


            disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(inp_folder_path, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(inp_folder_path, outputFileName{nn}, inputs_tblut,...
                inputs_tblut.RT.compute_reflectivity_uvSpec);


            % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
            idx_wl = source_wavelength>=(changing_variables(nn,end-2) - wl_perturb) &...
                source_wavelength<=(changing_variables(nn,end-1) + wl_perturb);

            [Refl_model_tblut(nn), ~] = reflectanceFunction_ver2(inputs_tblut, ds,...
                source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



        end


    else


        parfor nn = 1:num_INP_files
            % for ww = 1:size(inputs.RT.wavelengths2run, 1)


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(inp_folder_path, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(inp_folder_path, outputFileName{nn}, inputs_tblut,...
                inputs_tblut.RT.compute_reflectivity_uvSpec);


            % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
            idx_wl = source_wavelength>=(changing_variables(nn,end-2) - wl_perturb) &...
                source_wavelength<=(changing_variables(nn,end-1) + wl_perturb);

            [Refl_model_tblut(nn), ~] = reflectanceFunction_ver2(inputs_tblut, ds,...
                source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



        end


    end










    % save the calculated reflectances and the inputs
    save(folder_paths.saveOutput_filename, "inputs_tblut", "Refl_model_tblut"); % save inputSettings to the same folder as the input and output file



elseif inputs_tblut.flags.runUVSPEC == false

    load([inputs_tblut.savedCalculations_folderName,inputs_tblut.saveCalculations_fileName] ,'inputs_tblut','Refl_model_tblut');

end


%% ----- Find the minimum root-mean-square effective radius and optical depth -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

if isfield(simulated_measurements, 'Refl_model')==true

    tblut_retrieval = leastSquaresGridSearch_HySICS(simulated_measurements.Refl_model, Refl_model_tblut, inputs_tblut, folder_paths);

elseif isfield(simulated_measurements, 'Refl_model_2save')==true

    tblut_retrieval = leastSquaresGridSearch_HySICS(simulated_measurements.Refl_model_2save, Refl_model_tblut, inputs_tblut, folder_paths);

else

    error([newline, 'What simualted measurements should I use?', newline])

end



end
