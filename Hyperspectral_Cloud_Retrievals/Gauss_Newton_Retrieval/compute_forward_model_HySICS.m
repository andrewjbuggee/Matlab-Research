% ---- Compute fowrad model using MODIS spectral channels ----

% this function will compute reflectances using DISORT radiative
% transfer solver. It will do this for a water cloud with a droplet
% profile over any number spectral channels


% For the retrieval of r_top, r_bot, tau_c

% By Andrew J. Buggee
%%
function measurement_estimate = compute_forward_model_HySICS(current_guess, GN_inputs, spec_response, folder_paths)

disp([newline, 'Estimating spectral measurements...', newline])

% --- compute the forward model at our current estimate ---
r_top = current_guess(1);
r_bottom = current_guess(2);
tau_c = current_guess(3);



% ----- unpack parallel for loop variables ------
% We want to avoid large broadcast variables!
wavelengths2run = GN_inputs.RT.wavelengths2run;
libRadtran_inp = folder_paths.libRadtran_inp;
libRadtran_data_path = folder_paths.libRadtran_data;
wc_folder_path = folder_paths.libRadtran_water_cloud_files;
mie_folder_path = folder_paths.libRadtran_mie_folder;
which_computer = GN_inputs.which_computer;



% Read the solar flux file over the wavelength range specified
wavelength_vec = [min(wavelengths2run,[],"all"), max(wavelengths2run, [], "all")];

[source_flux, source_wavelength] = read_solar_flux_file(wavelength_vec, GN_inputs.RT.source_file);   % W/nm/m^2

% we will add and subtract a small fraction of the source file resolution
% to ensure rounding errors don't cause an issue when selecting the
% wavelengths needed from the source file
wl_perturb = GN_inputs.RT.source_file_resolution/3;   % nm

% ----------------------------------------------------------

% --------------------------------------------
% create water cloud file with droplet profile
% --------------------------------------------



% constraint - the physical constraint (string) - there are four
%       different string options for a physical constraint:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.



re = create_droplet_profile2([r_top, r_bottom], GN_inputs.RT.z, GN_inputs.RT.indVar, GN_inputs.RT.profile_type);     % microns - effective radius vector



% -----------------------------------
% ---- Write a Water Cloud file! ----
% -----------------------------------

% ------------------------------------------------------
% --------------------VERY IMPORTANT ------------------
% ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
% ------------------------------------------------------
loop_var = 0;

wc_filename = write_wc_file(re, tau_c, GN_inputs.RT.z_topBottom, GN_inputs.RT.lambda_forTau,...
    GN_inputs.RT.distribution_str, GN_inputs.RT.distribution_var, GN_inputs.RT.vert_homogeneous_str,...
    GN_inputs.RT.parameterization_str, GN_inputs.RT.indVar, false, GN_inputs.which_computer,...
    loop_var, 2, wc_folder_path, mie_folder_path);
wc_filename = wc_filename{1};

% ------------------------------------------------------
% ------------------------------------------------------
measurement_estimate = zeros(size(GN_inputs.RT.wavelengths2run,1), 1);

parfor ww = 1:size(wavelengths2run,1)
% for ww = 1:size(GN_inputs.RT.wavelengths2run,1)

    % define the input file name
    inputFileName = [num2str(mean(wavelengths2run(ww,:))), '_','nm_rTop_', num2str(r_top),...
        '_rBot_', num2str(r_bottom),'_tauC_', num2str(tau_c), '.INP'];

    outputFileName = ['OUTPUT_',inputFileName(1:end-4)];


    % ----- Write an INP file --------
    write_INP_file(libRadtran_inp, libRadtran_data_path, wc_folder_path, inputFileName, GN_inputs,...
        wavelengths2run(ww,:), wc_filename);


    % ----------------------------------------------------
    % --------------- RUN RADIATIVE TRANSFER -------------
    % ----------------------------------------------------


    % compute INP file
    runUVSPEC_ver2(libRadtran_inp, inputFileName, outputFileName,...
        which_computer);


    % read .OUT file
    % radiance is in units of mW/nm/m^2/sr
    [ds,~,~] = readUVSPEC_ver2(libRadtran_inp, outputFileName, GN_inputs,...
        GN_inputs.RT.compute_reflectivity_uvSpec);


    % compute the reflectance **NEED SPECTRAL RESPONSE INDEX***
    idx_wl = source_wavelength>=(wavelengths2run(ww,1) - wl_perturb) &...
        source_wavelength<=(wavelengths2run(ww,2) + wl_perturb);

    [measurement_estimate(ww), ~] = reflectanceFunction_ver2(GN_inputs, ds,...
        source_flux(idx_wl), spec_response(ww,:)');

end






end