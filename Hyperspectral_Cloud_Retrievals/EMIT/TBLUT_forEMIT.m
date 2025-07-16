%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using EMIT data


% ---------------- INPUTS ---------------
% ---------------------------------------
% (1) emit - EMIT data structure

% (2) spec_response - EMIT spectral response for each channel

% (3) emitDataFolder - this is the data folder where the EMIT data is
% located

% (4) folderpaths - this is a structure with all the folder paths needed to
% read ad store files





% By Andrew John Buggee

%%

function tblut_retrieval = TBLUT_forEMIT(emit, spec_response, emitDataFolder, folder_paths)


disp([newline, 'Computing the TBLUT retrieval...', newline])


%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs_tblut = create_emit_inputs_TBLUT(emitDataFolder, folder_paths, emit, spec_response);


%% SPECIAL CASES

% Values for 27_Jan_2024 - ** pixel [1242, 973] **
% override the cloud top height
inputs_tblut.RT.z_topBottom = [0.700, 0.500];    % km

% Values for 27_Jan_2024 - ** pixel [1242, 973] **
inputs_tblut.RT.modify_total_columnWaterVapor = true;             % don't modify the full column
inputs_tblut.RT.waterVapor_column = 36;   % mm - milimeters of water condensed in a column



%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs_tblut = define_source_for_EMIT(inputs_tblut, emit);




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

    % only jump on indexes where there is a unique r and tau pair

    parfor nn = 1:num_rEff*num_tauC
        % for nn = 1:num_rEff*num_tauC

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------


        temp = write_wc_file(changing_variables(2*nn -1, 1), changing_variables(2*nn -1,2),...
            inputs_tblut.RT.z_topBottom,inputs_tblut.RT.lambda_forTau, inputs_tblut.RT.distribution_str,...
            inputs_tblut.RT.distribution_var,inputs_tblut.RT.vert_homogeneous_str, inputs_tblut.RT.parameterization_str,...
            inputs_tblut.RT.indVar, inputs_tblut.compute_weighting_functions, inputs_tblut.which_computer, nn+(nn-1), 1);

        temp_names{nn} = temp{1};

    end

    % set the odd and even values to have the same file names
    wc_filename(1:num_wl:num_INP_files, 1) = temp_names;
    wc_filename(2:num_wl:num_INP_files, 1) = temp_names;


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


        inputFileName{nn} = [num2str(mean(wavelengths)), '_',...
            'nm_rEff_', num2str(changing_variables(nn,1)), '_tauC_', num2str(changing_variables(nn,2)), '_',...
            inputs_tblut.RT.atm_file(1:end-4),'.INP'];



        outputFileName{nn} = ['OUTPUT_',inputFileName{nn}(1:end-4)];


        % ------------------ Write the INP File --------------------
        write_INP_file(folder_paths.libRadtran_inp, folder_paths.libRadtran_data_path, inputFileName{nn}, inputs_tblut,...
            wavelengths, wc_filename{nn});


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


    % grab the values of the spectral response
    spec_response_value = spec_response.value;

    % store the reflectances
    Refl_model_tblut = zeros(num_INP_files, 1);


    parfor nn = 1:num_INP_files
        % for nn = 1:num_INP_files


        disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


        % ----------------------------------------------------
        % --------------- RUN RADIATIVE TRANSFER -------------
        % ----------------------------------------------------



        % compute INP file
        runUVSPEC_ver2(folder_paths.libRadtran_inp, inputFileName{nn}, outputFileName{nn},...
            inputs_tblut.which_computer);


        % read .OUT file
        % radiance is in units of mW/nm/m^2/sr
        [ds,~,~] = readUVSPEC_ver2(folder_paths.libRadtran_inp, outputFileName{nn}, inputs_tblut,...
            inputs_tblut.RT.compute_reflectivity_uvSpec);


        % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
        idx_wl = source_wavelength>=(changing_variables(nn,3) - wl_perturb) &...
            source_wavelength<=(changing_variables(nn,4) + wl_perturb);

        [Refl_model_tblut(nn), ~] = reflectanceFunction_ver2(inputs_tblut, ds,...
            source_flux(idx_wl), spec_response_value(changing_variables(nn,end),:)');


        %         [Refl_model_tblut(nn),~, inputs_tblut] = runReflectanceFunction_4EMIT(inputs_tblut,...
        %             names, emit.spec_response.value);



    end



    % save the calculated reflectances and the inputs
    save(inputs_tblut.saveOutput_fileName, "inputs_tblut", "Refl_model_tblut"); % save inputSettings to the same folder as the input and output file



elseif inputs_tblut.flags.runUVSPEC == false

    load([inputs_tblut.savedCalculations_folderName,inputs_tblut.saveCalculations_fileName] ,'inputs_tblut','Refl_model_tblut');

end







%% ----- Find the minimum root-mean-square effective radius and optical depth -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

tblut_retrieval = leastSquaresGridSearch_EMIT(emit.reflectance, Refl_model_tblut, inputs_tblut);



end
