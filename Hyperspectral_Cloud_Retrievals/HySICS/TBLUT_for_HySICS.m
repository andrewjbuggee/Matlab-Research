%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using simulated HySICS data


% ---------------- INPUTS ---------------
% ---------------------------------------
% (1) simulated_reflectance - HySICS simulated measurements

% (2) folderpaths - this is a structure with all the folder paths needed to
% read and store files, including the path to the simulated data set, the
% location to save the retirevals, and the location to save the INP files




% By Andrew John Buggee

%%

function tblut_retrieval = TBLUT_for_HySICS(simulated_reflectance, folder_paths)



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




   




idx = 0;


if inputs.flags.writeINPfiles == true



    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    % length of each independent variable
    num_rEff = length(inputs.RT.re);
    num_tauC = length(inputs.RT.tau_c);
    num_wl = length(inputs.bands2run);

    inputFileName = cell(num_rEff, num_tauC, 2);
    outputFileName = cell(num_rEff, num_tauC, 2);


    for rr = 1:num_rEff



        for tc = 1:num_tauC


            idx = idx + 1;
            % -----------------------------------
            % ---- Write a Water Cloud file! ----
            % -----------------------------------

            % ------------------------------------------------------
            % --------------------VERY IMPORTANT ------------------
            % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
            % ------------------------------------------------------
            wc_filename = write_wc_file(inputs.RT.re(rr), inputs.RT.tau_c(tc), inputs.RT.z_topBottom,...
                inputs.RT.lambda_forTau, inputs.RT.distribution_str, inputs.RT.distribution_var,...
                inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str,...
                inputs.which_computer, idx);

            inputs.RT.wc_filename = wc_filename{1};




            % step through each band
            for ww = 1:num_wl


                % set the wavelengths for each file
                inputs.RT.wavelengths = inputs.RT.wavelengths2run(ww,:);
                

                % ------------------------------------------------
                % ---- Define the input and output filenames! ----
                % ------------------------------------------------
                % input_names need a unique identifier. Let's give them the nn value so
                % they can be traced, and are writen over in memory


                inputFileName{rr, tc, ww} = [num2str(mean(inputs.RT.wavelengths)), '_',...
                    'nm_rEff_', num2str(inputs.RT.re(rr)), '_tauC_', num2str(inputs.RT.tau_c(tc)), '_',...
                    inputs.RT.atm_file(1:end-4),'.INP'];



                outputFileName{rr, tc, ww} = ['OUTPUT_',inputFileName{rr,tc,ww}(1:end-4)];


                % ------------------ Write the INP File --------------------
                write_INP_file(folder_paths.libRadtran_inp, inputs.libRadtran_data_path, inputFileName{rr, tc, ww}, inputs);




            end



        end

    end






else

    % if the files already exist, just grab the names!
    [names.inp, inputs] = getMODIS_INPnames_withClouds(simulated_reflectance.solar, inputs, pixels2use);
    names.out = writeOutputNames(names.inp);

end



%%

% store the reflectances
Refl_model = zeros(num_wladfdsf, num_rEff, num_tauC);


for rr = 1:num_rEff



    for tc = 1:num_tauC





        parfor ww = 1:num_wl
            % for ww = 1:size(inputs.RT.wavelengths2run, 1)


            disp(['Iteration: [re, tc, ww] = [', num2str(rr), '/', num2str(num_rEff),', ',...
                num2str(tc), '/', num2str(num_tauC), ', ', num2str(ww),'/',num2str(num_wl),...
                ']', newline])


            % ----------------------------------------------------
            % --------------- RUN RADIATIVE TRANSFER -------------
            % ----------------------------------------------------


            % compute INP file
            [inputSettings] = runUVSPEC(inputs.folderpath_inp, inputFileName{rr, tc, ww}, outputFileName{rr, tc, ww});

            % read .OUT file
            % radiance is in units of mW/nm/m^2/sr
            [ds,~,~] = readUVSPEC(inputs.folderpath_inp, outputFileName{rr, tc, ww},inputSettings(2,:),...
                inputs.RT.compute_reflectivity_uvSpec);

            % Store the Radiance
            %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

            % compute the reflectance
            % Refl_model(ww, rr, tc) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
            %     spec_response.value(ww, :)');

            [Refl_model(ww, rr, tc), ~] = reflectanceFunction(inputSettings(2,:), ds, spec_response.value(ww,:));








        end



    end

end



%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs.flags.runUVSPEC == true

    % 1st output - R is the reflectance integrated over a bandwidth
    % 2nd output - Rl is the reflectance at each spectral bin
    tic
    [R,~, inputs] = runReflectanceFunction_4EMIT(inputs, names, simulated_reflectance.spec_response.value);
    toc

    % Save the pixels2use structure
    save([inputs.folder2save.reflectance_calcs, inputs.reflectance_calculations_fileName],...
        "pixels2use", "-append"); % save inputSettings to the same folder as the input and output file

elseif inputs.flags.runUVSPEC == false

    load([inputs.savedCalculations_folderName,inputs.saveCalculations_fileName] ,'inputs','R');

end


%% ----- Find the minimum root-mean-square effective radius and optical depth -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

tblut_retrieval = leastSquaresGridSearch_EMIT(simulated_reflectance.reflectance, R, inputs);



end
