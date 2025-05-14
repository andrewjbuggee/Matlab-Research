%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using simulated HySICS data


% ---------------- INPUTS ---------------
% ---------------------------------------
% (1) simulated_reflectance - HySICS simulated measurements

% (2) folderpaths - this is a structure with all the folder paths needed to
% read and store files, including the path to the simulated data set, the
% location to save the retirevals, and the location to save the INP files




% By Andrew John Buggee

%%

function tblut_retrieval = TBLUT_for_HySICS_ver2(simulated_reflectance, folder_paths)



%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_HySICS_inputs_TBLUT(folder_paths);



%% Find the measurements closest to the bands to run

[~, idx_1] = min(abs(simulated_reflectance.inputs.bands2run - inputs.bands2run(1)));

[~, idx_2] = min(abs(simulated_reflectance.inputs.bands2run - inputs.bands2run(2)));

% error if the values found are at least 15nm from the intended wavlengths
if abs(mean(simulated_reflectance.inputs.RT.wavelengths2run(idx_1,:)) - 650)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 650 nm', newline])

elseif abs(mean(simulated_reflectance.inputs.RT.wavelengths2run(idx_2,:)) - 2130)>15

    error([newline, 'The measurements provided dont have a reflectance measurement close to 650 nm', newline])

else

    % Then we set the bands to run to be to ones found to be closest to the
    % desired bands out of the measurement bands provided
    inputs.bands2run_from_set_of_measurements = [idx_1, idx_2];
    inputs.bands2plot = inputs.bands2run;

    % ---- Define the wavelengths ----
    inputs.RT.wavelengths2run = simulated_reflectance.inputs.RT.wavelengths2run(inputs.bands2run_from_set_of_measurements,:);



end
%% ----- Create .INP files for HySICS TBLUT -----


if inputs.flags.writeINPfiles == true



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

    parfor nn = 1:num_rEff*num_tauC

        % -----------------------------------
        % ---- Write a Water Cloud file! ----
        % -----------------------------------


        temp = write_wc_file(changing_variables(2*nn -1, 1), changing_variables(2*nn -1,2),...
            inputs.RT.z_topBottom,inputs.RT.lambda_forTau, inputs.RT.distribution_str,...
            inputs.RT.distribution_var,inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
            inputs.which_computer, nn+(nn-1));

        temp_names{nn} = temp{1};

    end

    % set the odd and even values to have the same file names
    wc_filename(1:num_wl:num_INP_files, 1) = temp_names;
    wc_filename(2:num_wl:num_INP_files, 1) = temp_names;


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
        write_INP_file(folder_paths.libRadtran_inp, inputs.libRadtran_data_path, inputFileName{nn}, inputs,...
            wavelengths, wc_filename{nn});


    end




else

    % if the files already exist, just grab the names!
    error([newline, 'Dont know how to do this!', newline])

    % [names.inp, inputs] = getMODIS_INPnames_withClouds(simulated_reflectance.solar, inputs, pixels2use);
    % names.out = writeOutputNames(names.inp);

end







%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs.flags.runUVSPEC == true

    spec_response = simulated_reflectance.spec_response.value;


    % store the reflectances
    Refl_model = zeros(num_INP_files, 1);


    parfor nn = 1:num_INP_files
        % for ww = 1:size(inputs.RT.wavelengths2run, 1)


        disp(['Iteration: nn/total_files = [', num2str(nn), '/', num2str(num_INP_files),']', newline])


        % ----------------------------------------------------
        % --------------- RUN RADIATIVE TRANSFER -------------
        % ----------------------------------------------------


        % compute INP file
        [inputSettings] = runUVSPEC(folder_paths.libRadtran_inp, inputFileName{nn}, outputFileName{nn},...
            inputs.which_computer);

        % read .OUT file
        % radiance is in units of mW/nm/m^2/sr
        [ds,~,~] = readUVSPEC(folder_paths.libRadtran_inp, outputFileName{nn},inputSettings(2,:),...
            inputs.RT.compute_reflectivity_uvSpec);

        % Store the Radiance
        %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

        % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
        [Refl_model(nn), ~] = reflectanceFunction(inputSettings(2,:), ds,...
            spec_response(changing_variables(nn,end),:));



    end



    % save the calculated reflectances and the inputs
    save(inputs.save_mat_filename, "inputs", "Refl_model"); % save inputSettings to the same folder as the input and output file



elseif inputs.flags.runUVSPEC == false

    load([inputs.savedCalculations_folderName,inputs.saveCalculations_fileName] ,'inputs','Refl_model');

end


%% ----- Find the minimum root-mean-square effective radius and optical depth -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

tblut_retrieval = leastSquaresGridSearch_HySICS(simulated_reflectance.Refl_model, Refl_model, inputs);



end
