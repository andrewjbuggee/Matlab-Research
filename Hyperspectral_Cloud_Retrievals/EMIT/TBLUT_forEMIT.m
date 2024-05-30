%% Compute the Two-Wavelength Estimate of Effective radius and optical depth using EMIT data


% ---------------- INPUTS ---------------
% ---------------------------------------
% (1) emit - EMIT data structure

% (2) emitDataFolder - this is the data folder where the EMIT data is
% located

% (3) pixels2use - these are the pixels to use for the retrieval



% By Andrew John Buggee

%%

function tblut_retrieval = TBLUT_forEMIT(emit, emitDataFolder, folder2save, pixels2use)



%% Create an input structure that helps write the INP files

% this is a built-in function that is defined at the bottom of this script
inputs = create_emit_inputs_TBLUT(emitDataFolder, folder2save, emit);


%% Define the solar source file name and read in the solar source data

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


%% ----- Create .INP files for EMIT TBLUT -----


if inputs.flags.writeINPfiles == true

    [names.inp, inputs] = write_INP_file_4EMIT_homogenous(inputs, pixels2use, emit);

    % now lets write the output names

    names.out = writeOutputNames(names.inp);
else

    % if the files already exist, just grab the names!
    [names.inp, inputs] = getMODIS_INPnames_withClouds(emit.solar, inputs, pixels2use);
    names.out = writeOutputNames(names.inp);

end


%% ----- Run uvspec and calculate Reflectance Function using LibRadTran -----

% geometry stays the same, but we calculate the radiative transfer equation
% for different values of effective radius and optical depth

if inputs.flags.runUVSPEC == true

    % 1st output - R is the reflectance integrated over a bandwidth
    % 2nd output - Rl is the reflectance at each spectral bin
    tic
    [R,~] = runReflectanceFunction_4EMIT(inputs, names, emit.spec_response);
    toc

elseif inputs.flags.runUVSPEC == false

    load([inputs.savedCalculations_folderName,inputs.saveCalculations_fileName] ,'inputs','R');

end


%% ----- Find the minimum root-mean-square effective radius and optical depth -----

% first grid search is on a coarse grid
% we want to minimize two the reflectance for two wavelengths

% if interpGridScalFactor is 10, then 9 rows will be interpolated to be 90
% rows, and 10 columns will be interpolated to be 100 columns

tblut_retrieval = leastSquaresGridSearch_EMIT(emit.reflectance, R, inputs);



end
