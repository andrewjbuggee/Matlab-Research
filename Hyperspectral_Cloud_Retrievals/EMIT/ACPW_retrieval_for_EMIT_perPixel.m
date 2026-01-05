%% Estimate the above cloud precipitable water amount using Simualted HySICS measurements


% By Andrew John Buggee

%%

function [acpw_retrieval] = ACPW_retrieval_for_EMIT_perPixel(emit, spec_response, tblut_retrieval, folder_paths, use_MODIS_AIRS_data,...
    GN_inputs, print_libRadtran_err, print_status_updates, pixel_num)



%% Unpack folder_paths

% define the INP folder location
libRadtran_inp = folder_paths.libRadtran_inp;

% define the libRadtran data path
libRadtran_data_path = folder_paths.libRadtran_data;

% define the water cloud directory
libRadtran_water_cloud_files = folder_paths.libRadtran_water_cloud_files;

% define the mie calculations directory
libRadtran_mie_folder = folder_paths.libRadtran_mie_folder;

% define the atmmod folder where the custom profiles are stored
atm_folder_path = folder_paths.atm_folder_path;

which_computer = folder_paths.which_computer;


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_acpw = create_emit_inputs_ACPW_perPixel(emit, tblut_retrieval, print_libRadtran_err, pixel_num);


%% Update based on GN_inputs

if exist("use_MODIS_AIRS_data", "var")==1 && use_MODIS_AIRS_data==true

    % override the cloud top height
    inputs_acpw.RT.z_topBottom = GN_inputs.RT.z_topBottom;  % km

    % change inputs that depend on z_topBottom
    % Water Cloud depth
    inputs_acpw.RT.H = inputs_acpw.RT.z_topBottom(1) - inputs_acpw.RT.z_topBottom(2);                                % km - geometric thickness of cloud


    inputs_acpw.RT.use_radiosonde_file = true;

    % For the ACPW retrieval, we will use the just temperature and
    % pressure from the radiosonde file
    % ** only use temp and pres **
    inputs_acpw.RT.radiosonde_file = GN_inputs.RT.radiosonde_file_T_P;

end

%% Find the measurements closest to the bands to run

if size(emit.radiance.measurements, 1) == 285

    bands = 1:285;

[~, idx_1] = min(abs(bands - inputs_acpw.bands2run(1)));

[~, idx_2] = min(abs(bands - inputs_acpw.bands2run(2)));

[~, idx_3] = min(abs(bands - inputs_acpw.bands2run(3)));

else

    error([newline, 'Have you downselected the number of EMIT wavelengths to use?', newline])

end


% error if the values found are at least 15nm from the intended wavlengths
if abs(emit.radiance.wavelength(idx_1) - 900)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 900 nm', newline])

elseif abs(emit.radiance.wavelength(idx_2) - 955)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 955 nm', newline])


elseif abs(emit.radiance.wavelength(idx_3) - 1127)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 1127 nm', newline])

else

    % Then we set the bands-to-run to be to ones found to be closest to the
    % desired bands out of the measurement bands provided
    inputs_acpw.bands2run_from_set_of_measurements = [idx_1, idx_2, idx_3];
    inputs_acpw.bands2plot = inputs_acpw.bands2run;

    % ---- Define the wavelengths ----
    inputs_acpw.RT.wavelengths2run = [spec_response.wavelength(inputs_acpw.bands2run_from_set_of_measurements, 1),...
        spec_response.wavelength(inputs_acpw.bands2run_from_set_of_measurements, end)];



end





%% ----- Create .INP files for HySICS ACPW -----


if inputs_acpw.flags.writeINPfiles == true



    inputs_acpw.acpw_sim = 0:40;    % mm

    % num wavelengths
    num_wl = length(inputs_acpw.bands2run);

    % length of each independent variable
    num_tcpw = length(inputs_acpw.acpw_sim);


    num_INP_files = num_tcpw * num_wl;

    inputFileName = cell(num_INP_files, 1);
    outputFileName = cell(num_INP_files, 1);


    % changing variable steps through tcpw and wavelength
    % in for loop speak, it would be:
    % for pw = 1:num_tcpw
    %   for ww = 1:num_wl
    changing_variables = [reshape(repmat(inputs_acpw.acpw_sim, num_wl,1), [],1),...
        repmat(inputs_acpw.RT.wavelengths2run, num_tcpw, 1)];


    % Add a final column that includes the index for the spectral response
    % function. These always increase chronologically
    changing_variables = [changing_variables, repmat((1:3)',  num_tcpw, 1)];
    % changing_variables = [changing_variables, repmat(inputs_acpw.bands2run_from_set_of_measurements',  num_tcpw, 1)];

    % Write the water cloud file


    wc_filename = write_wc_file(tblut_retrieval.minRe, tblut_retrieval.minTau,...
        inputs_acpw.RT.z_topBottom,inputs_acpw.RT.lambda_forTau, inputs_acpw.RT.distribution_str,...
        inputs_acpw.RT.distribution_var, inputs_acpw.RT.vert_homogeneous_str, inputs_acpw.RT.parameterization_str,...
        inputs_acpw.RT.indVar, inputs_acpw.compute_weighting_functions, which_computer,...
        1, inputs_acpw.RT.num_re_parameters, libRadtran_water_cloud_files,...
        libRadtran_mie_folder);

    wc_filename = wc_filename{1};






    % Now write all the INP files
    parfor nn = 1:num_INP_files
    % for nn = 1:num_INP_files


        % set the wavelengths for each file
        wavelengths = changing_variables(nn, end-2:end-1);

        % create a custom water vapor profile
        custom_waterVapor_profile = alter_aboveCloud_columnWaterVapor_profile(inputs_acpw, changing_variables(nn,1),...
            atm_folder_path);

        % ------------------------------------------------
        % ---- Define the input and output filenames! ----
        % ------------------------------------------------
        % input_names need a unique identifier. Let's give them the nn value so
        % they can be traced, and are writen over in memory


        inputFileName{nn} = ['ACPW_retrieval_', num2str(mean(wavelengths)), '_','nm_re_', num2str(round(tblut_retrieval.minRe, 3)),...
            '_tauC_', num2str(round(tblut_retrieval.minTau, 3)), '_acpw_',...
            num2str(round(changing_variables(nn,1), 3)),'mm.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(libRadtran_inp, libRadtran_data_path, libRadtran_water_cloud_files,...
            inputFileName{nn}, inputs_acpw, wavelengths, wc_filename, [], tblut_retrieval.minTau,...
            custom_waterVapor_profile, []);


    end



else

    % if the files already exist, just grab the names!
    error([newline, 'Dont know how to do this!', newline])


end


%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs_acpw.flags.runUVSPEC == true


    % Read the solar flux file over the wavelength range specified
    wavelength_vec = [min(inputs_acpw.RT.wavelengths2run,[],"all"), max(inputs_acpw.RT.wavelengths2run, [], "all")];

    [source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, inputs_acpw.RT.source_file);   % W/nm/m^2

    % we will add and subtract a small fraction of the source file resolution
    % to ensure rounding errors don't cause an issue when selecting the
    % wavelengths needed from the source file
    wl_perturb = inputs_acpw.RT.source_file_resolution/3;   % nm




    % Let's only keep values we need
    spec_response = spec_response.value(inputs_acpw.bands2run_from_set_of_measurements, :);


    % store the reflectances
    Refl_model_acpw = zeros(num_INP_files, 1);


    if print_status_updates==true

        parfor nn = 1:num_INP_files
            % for ww = 1:size(inputs.RT.wavelengths2run, 1)


            disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs_acpw,...
                inputs_acpw.RT.compute_reflectivity_uvSpec);


            % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
            idx_wl = source_wavelength>=(changing_variables(nn,end-2) - wl_perturb) &...
                source_wavelength<=(changing_variables(nn,end-1) + wl_perturb);

            [Refl_model_acpw(nn), ~] = reflectanceFunction_ver2(inputs_acpw, ds,...
                source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



        end


    else


        parfor nn = 1:num_INP_files
            % for ww = 1:size(inputs.RT.wavelengths2run, 1)


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------



            % compute INP file
            runUVSPEC_ver2(libRadtran_inp, inputFileName{nn}, outputFileName{nn},...
                which_computer);


            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName{nn}, inputs_acpw,...
                inputs_acpw.RT.compute_reflectivity_uvSpec);


            % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
            idx_wl = source_wavelength>=(changing_variables(nn,end-2) - wl_perturb) &...
                source_wavelength<=(changing_variables(nn,end-1) + wl_perturb);

            [Refl_model_acpw(nn), ~] = reflectanceFunction_ver2(inputs_acpw, ds,...
                source_flux(idx_wl), spec_response(changing_variables(nn,end),:)');



        end


    end




elseif inputs_acpw.flags.runUVSPEC == false

    load([inputs_acpw.savedCalculations_folderName,inputs_acpw.saveCalculations_fileName] ,'inputs_acpw','Refl_model_acpw');

end




%% Find the minimum RMS difference between the measurements and the calculations



R_measurement = emit.reflectance.value(inputs_acpw.bands2run_from_set_of_measurements, pixel_num);

RMS = sqrt( mean( (repmat(R_measurement, 1, num_tcpw) - reshape(Refl_model_acpw, num_wl, [])).^2, 1) );

[~, idx_min] = min(RMS);

acpw_retrieval.min_nonInterpolated = inputs_acpw.acpw_sim(idx_min);


%% Calculations show that reflectance decreases monotonically with increasing above cloud column water vapor
% Therefore, we can linearly inertoplate to get a more accurate number

inputs_acpw.acpw_interp1 = (min(inputs_acpw.acpw_sim):0.01:max(inputs_acpw.acpw_sim));

Refl_model_acpw_interp1 = [interp1(inputs_acpw.acpw_sim, Refl_model_acpw(1:3:end)', inputs_acpw.acpw_interp1);...
                         interp1(inputs_acpw.acpw_sim, Refl_model_acpw(2:3:end)', inputs_acpw.acpw_interp1);...
                         interp1(inputs_acpw.acpw_sim, Refl_model_acpw(3:3:end)', inputs_acpw.acpw_interp1)];

% Find the new minimum RMS

RMS_interp = sqrt( mean( (repmat(R_measurement, 1, numel(inputs_acpw.acpw_interp1)) - Refl_model_acpw_interp1).^2, 1) );

[~, idx_min_interp] = min(RMS_interp);

acpw_retrieval.min_interpolated = inputs_acpw.acpw_interp1(idx_min_interp);



%% save the calculated reflectances and the inputs
if isfile(folder_paths.saveOutput_filename)==true

    save(folder_paths.saveOutput_filename, "inputs_acpw", "Refl_model_acpw", "acpw_retrieval", '-append'); % save inputSettings to the same folder as the input and output file

else

    save(folder_paths.saveOutput_filename, "inputs_acpw", "Refl_model_acpw", "acpw_retrieval"); % save inputSettings to the same folder as the input and output file

end

end