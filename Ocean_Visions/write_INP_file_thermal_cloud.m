%% Write a libRadtran (uvspec) input file for a broadband thermal (longwave) calculation
%
% Builds an .INP file that computes downwelling and upwelling longwave
% irradiance at a specified altitude for a clear sky, a single-layer liquid
% cloud, or a single-layer ice cloud.
%
% Design notes -- why these options:
%
%   source thermal          Emission-only calculation. There is no solar term at
%                           all, which is the right idealization for the Arctic in
%                           late December (polar night poleward of the Arctic
%                           circle on 22 Dec; the solar zenith angle never drops
%                           below 90 deg, so the shortwave source is exactly zero).
%                           Note that day_of_year / time only enter uvspec through
%                           the Earth-Sun distance for the *solar* source, so they
%                           are deliberately omitted -- the December condition is
%                           carried entirely by the choice of atmosphere_file.
%
%   mol_abs_param reptran   Representative-wavelength band parameterization
%                           (Gasteiger et al., 2014). The thermal set spans
%                           2501.6 - 99808.6 nm. 'coarse' is the standard broadband
%                           flux resolution; 'medium'/'fine' are also available.
%
%   output_process per_band Each output row is a band-integrated irradiance in
%                           W/m^2 (the thermal default unit is W/(m^2 cm^-1), and
%                           per_band converts to W/m^2 per band). Summing the
%                           column in MATLAB then gives the broadband longwave
%                           irradiance AND keeps the spectrum for diagnostics --
%                           which output_process sum would throw away.
%
%   wc_properties hu        Hu and Stamnes (1993). Valid 0.29 - 150 micron, so it
%                           covers the full thermal range. It is a flux-oriented
%                           parameterization (Henyey-Greenstein phase function),
%                           which is exactly the intended use here.
%
%   ic_properties yang      Key et al. (2002) / Yang; 0.2 - 100 micron, habit
%                           selectable via ic_habit. 'solid-column' is valid for
%                           r_eff in [5.96, 84.22] micron.
%                           'fu' (Fu 1996; Fu et al. 1998) is the libRadtran
%                           default and also assumes randomly oriented hexagonal
%                           columns; it agrees with yang/solid-column to better
%                           than 0.1 W/m^2 for these cases, and is a useful
%                           cross-check.
%
%   ic_properties yang2013  Yang et al. (2013); 0.2 - 99 micron, 9 habits x 3
%                           roughness levels, r_eff 5 - 90 micron for every
%                           habit. Selected with ic_habit_yang2013 (habit AND
%                           roughness). It can only be used together with
%                           inputs.wavelength_index -- see the spectral range
%                           block below for why. baum_v36 (tables also stop at
%                           99000 nm) would need the same treatment.
%
%   albedo                  In the thermal, uvspec sets the surface emissivity to
%                           (1 - albedo). Surface properties affect the upwelling
%                           flux directly and the downwelling flux only through
%                           the small reflected-then-rescattered term.
%
%
% INPUTS:
%   (1) INP_folderpath - folder the .INP file is written into. Created if absent.
%
%   (2) inputFileName - name of the .INP file (including the extension).
%
%   (3) inputs - structure holding the run settings. Required fields:
%
%       inputs.libRadtran_data_path - full path to libRadtran's data/ folder
%       inputs.atmosphere_file      - full path to the atmospheric profile
%       inputs.rte_solver           - e.g. 'twostr'
%       inputs.mol_abs_param        - e.g. 'reptran coarse'
%       inputs.wavelength_bounds_nm - [lambda_min_nm, lambda_max_nm]
%       inputs.albedo               - surface albedo (thermal emissivity = 1 - albedo)
%       inputs.zout_km              - output altitude, km above the surface
%       inputs.cloud_phase          - 'none', 'liquid', or 'ice'
%       inputs.cloud_file_path      - full path to the wc/ic .DAT file
%                                     (ignored when cloud_phase is 'none')
%       inputs.wc_properties        - e.g. 'hu'          (liquid clouds only)
%       inputs.ic_properties        - e.g. 'yang'        (ice clouds only)
%       inputs.ic_habit             - e.g. 'solid-column' (ice clouds only; set
%                                     to [] or '' to fall back to the default habit)
%
%       Optional fields:
%
%       inputs.wavelength_index     - [idx_lo, idx_hi] band indices. When present
%                                     and non-empty this REPLACES the wavelength
%                                     option (uvspec rejects both together).
%       inputs.ic_habit_yang2013    - e.g. 'solid_column'. REQUIRED when
%                                     ic_properties is yang2013.
%       inputs.ic_roughness         - 'smooth', 'moderate', or 'severe'. REQUIRED
%                                     when ic_properties is yang2013.
%
%
% OUTPUTS:
%   (1) An .INP file written to INP_folderpath.
%
%
% By Andrew John Buggee
%%

function [] = write_INP_file_thermal_cloud(INP_folderpath, inputFileName, inputs)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

if nargin~=3
    error([newline, 'Need 3 inputs: the INP folder path, the INP file name, and the inputs structure.', newline])
end

required_fields = {'libRadtran_data_path', 'atmosphere_file', 'rte_solver', 'mol_abs_param',...
    'wavelength_bounds_nm', 'albedo', 'zout_km', 'cloud_phase'};

for nn = 1:length(required_fields)
    if isfield(inputs, required_fields{nn})==false
        error([newline, 'The inputs structure is missing the field: ', required_fields{nn}, newline])
    end
end

if exist(INP_folderpath, 'dir')==0
    mkdir(INP_folderpath)
end

if INP_folderpath(end)~=filesep
    INP_folderpath = [INP_folderpath, filesep];
end


% ------------------------------------------------------------
% -------------------- WRITE THE .INP FILE -------------------
% ------------------------------------------------------------

fileID = fopen([INP_folderpath, inputFileName], 'w');
if fileID<0
    error([newline, 'Could not open the INP file for writing: ', INP_folderpath, inputFileName, newline])
end


% --- where libRadtran finds its own data ---
fprintf(fileID, '%s %s \n', 'data_files_path', inputs.libRadtran_data_path);
fprintf(fileID, '%s %s \n\n', 'atmosphere_file', inputs.atmosphere_file);

% --- radiation source: thermal emission only, no solar term ---
fprintf(fileID, '%s \n', 'source thermal');

% --- gas absorption parameterization ---
fprintf(fileID, '%s %s \n', 'mol_abs_param', inputs.mol_abs_param);

% --- radiative transfer equation solver ---
fprintf(fileID, '%s %s \n\n', 'rte_solver', inputs.rte_solver);

% --- radiative transfer streams, if needed ---
if strcmp(inputs.rte_solver, 'disort')==true

    fprintf(fileID, '%s %s \n\n', 'number_of_streams', num2str(inputs.n_streams));
end


% --- spectral range ---
% Two mutually exclusive ways to select the spectral range. uvspec aborts with
% "it does not make sense to define both 'wavelength' and 'wavelength_index'"
% if both appear, so wavelength_index takes precedence when it is supplied.
%
% wavelength_index selects bands of the mol_abs_param parameterization by
% index rather than by wavelength. It is the only way to drop the topmost
% REPTRAN thermal band (index 260 of 260 for 'reptran coarse', centred at
% 93478 nm and extending to 99808.6 nm). That matters because uvspec validates
% an ice optical property table against the FULL band set of the chosen
% parameterization, not against the subset requested by 'wavelength' -- so
% ic_properties yang2013 (tables stop at 99000 nm) is rejected no matter what
% bounds 'wavelength' is given. Dropping band 260 makes yang2013 usable.
if isfield(inputs, 'wavelength_index')==true && isempty(inputs.wavelength_index)==false

    fprintf(fileID, '%s %d %d \n', 'wavelength_index', inputs.wavelength_index(1),...
        inputs.wavelength_index(2));

else

    fprintf(fileID, '%s %f %f \n', 'wavelength', inputs.wavelength_bounds_nm(1),...
        inputs.wavelength_bounds_nm(2));

end

% --- surface: albedo, so thermal emissivity = 1 - albedo ---
fprintf(fileID, '%s %f \n', 'albedo', inputs.albedo);

% --- output altitude, km above the surface ---
fprintf(fileID, '%s %f \n\n', 'zout', inputs.zout_km);


% --- the cloud ---
if strcmp(inputs.cloud_phase, 'liquid')==true

    fprintf(fileID, '%s %s \n', 'wc_file 1D', inputs.cloud_file_path);
    fprintf(fileID, '%s %s \n\n', 'wc_properties', inputs.wc_properties);

elseif strcmp(inputs.cloud_phase, 'ice')==true

    fprintf(fileID, '%s %s \n', 'ic_file 1D', inputs.cloud_file_path);
    fprintf(fileID, '%s %s \n', 'ic_properties', inputs.ic_properties);

    % Yang et al. (2013) uses its own habit keyword, which also requires a
    % surface roughness (smooth / moderate / severe). Every other habit-aware
    % parameterization (key, yang, hey, baum_v36) uses plain ic_habit.
    if contains(inputs.ic_properties, 'yang2013')==true

        if isfield(inputs, 'ic_habit_yang2013')==false || isempty(inputs.ic_habit_yang2013)==true
            error([newline, 'ic_properties yang2013 requires inputs.ic_habit_yang2013.', newline])
        end

        if isfield(inputs, 'ic_roughness')==false || isempty(inputs.ic_roughness)==true
            error([newline, 'ic_properties yang2013 requires inputs.ic_roughness',...
                ' (smooth, moderate, or severe).', newline])
        end

        fprintf(fileID, '%s %s %s \n', 'ic_habit_yang2013', inputs.ic_habit_yang2013,...
            inputs.ic_roughness);

    elseif isfield(inputs, 'ic_habit')==true && isempty(inputs.ic_habit)==false

        fprintf(fileID, '%s %s \n', 'ic_habit', inputs.ic_habit);

    end

    fprintf(fileID, '\n');

elseif strcmp(inputs.cloud_phase, 'none')==true

    % clear sky - no cloud lines written

else

    error([newline, 'inputs.cloud_phase must be ''none'', ''liquid'', or ''ice''.', newline])

end


% --- what to write out ---
% for the twostr solver the natural columns are lambda, edir, edn, eup, uavg.
% edir is identically zero for a thermal-only calculation, but it is kept so
% that the column layout is explicit and self-documenting.
fprintf(fileID, '%s \n', 'output_user lambda edir edn eup');

% band-integrated irradiance, W/m^2 per reptran band
fprintf(fileID, '%s \n\n', 'output_process per_band');

% do you want to print the error messages?
fprintf(fileID, '%s \n', inputs.err_msg);

fclose(fileID);


end
