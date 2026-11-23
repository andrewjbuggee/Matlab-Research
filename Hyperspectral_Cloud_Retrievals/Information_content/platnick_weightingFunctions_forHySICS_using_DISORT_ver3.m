%% Compute Platnick (2000) maximum-penetration weighting functions for HySICS using DISORT - ver3
%
% By Andrew John Buggee (ver3 written with Claude, Nov 2026, from ver2)
%
% -------------------------------------------------------------------------
% WHAT CHANGED RELATIVE TO ver2
% -------------------------------------------------------------------------
% ver3 computes R(tau) from TRUNCATED clouds: the top tau of the full tau_c
% cloud, with every layer below tau removed. That is how Platnick (2000)
% defines the maximum-penetration weighting:
%
%   Platnick (2000), JGR 105(D18), section 4.1, eq. (4):
%       w_m(tau, tau_c) = [dR(tau)/dtau] / R(tau_c)
%   with R(tau) "the reflectance from the portion of the cloud above the
%   level tau (i.e., calculated with lower layers absent)". He adds that
%   "the derivative is calculated for a dtau layer added to the base of the
%   cloud, not the top, an important distinction when cloud microphysics
%   varies with height."
%
% ver2 instead wrote, for every tau in tau_2run, the FULL r_bot -> r_top
% profile with its LWC rescaled so the whole cloud had optical depth tau.
% A cloud with tau = 0.5 therefore still contained 5 um droplets at its
% base, although the top 0.5 of the real tau_c = 8 cloud only holds
% droplets of ~11.8-12 um.
%
% -------------------------------------------------------------------------
% WHY THE JUNE 2025 TRUNCATION RUN GAVE NEGATIVE VALUES
% -------------------------------------------------------------------------
% disort_HySICS_reflectance_for_weightingFunctions_inhomogeneous_droplet_profile_sim-ran-on-09-Jun-2025_rev1.mat
% (platnick_weightingFunctions_forHySICS_using_DISORT.m, layer removal via
% write_wc_file) has negative w in 18% of its values, the worst at -7.7
% times its band's peak. In 1902 of 16434 cases R DROPPED when a layer was
% added at the base (median drop 1.7e-3, max 1e-2, i.e. ~1-6% of R), mostly
% for thin truncated clouds (the top ~15 of 250 layers). Over a black
% surface adding a layer at the base can only increase R, so these are
% numerical errors, not physics.
% That run used sza = 0 and vza = 0 (exact backscatter). The
% 66Bands_zOut_at_cloudTop_20SZA run from the same week, with the same
% truncation code and cloud but sza = 20, vza = 10, has NO decreases in R.
% So the problem is tied to the exact-backscatter geometry, not to the
% truncation. My inference (not verified): DISORT's single-scattering
% (Nakajima-Tanaka) correction evaluates the phase function at a scattering
% angle of 180 deg, inside the glory, where it is most sensitive to
% droplet size and hardest to reconstruct from a finite set of Legendre
% moments. Thin clouds, whose reflectance is mostly single scattering,
% would be hit hardest - which is where the errors are.
% => Avoid sza = vza = 0 (a warning is issued below).
% ver3 additionally uses layers of EQUAL optical thickness, so no layer
% is optically tiny (that run's uniform-z layers ranged from 0.0005 to
% 0.054 in optical thickness), and never uses "wc_modify tau set"
% (enforced below), which would rescale the LWC of every retained layer
% and undo the truncation.
%
% -------------------------------------------------------------------------
% HOW ver3 BUILDS THE CLOUD
% -------------------------------------------------------------------------
% 1) The physical cloud is the one ver2 uses (built by write_wc_file):
%    LWC increases linearly with height from zero at cloud base, r_e(z)
%    comes from create_droplet_profile2 (adiabatic by default), the
%    geometric extent is inputs.RT.z_topBottom, and the total optical depth
%    at inputs.RT.lambda_forTau is tau_c. Keeping the same cloud isolates
%    the effect of truncation.
% 2) The cloud is re-layered into N_layers homogeneous layers of EQUAL
%    optical thickness dtau_layer. Platnick (2000, section 5) also built his
%    clouds from homogeneous layers of equal optical thickness (0.25).
%    Equal-tau layers give the same tau resolution at cloud top, where the
%    absorbing-band weighting functions peak, as at cloud base. Uniform-z
%    layers (what write_wc_file makes) are optically THICKEST at cloud top
%    because LWC is largest there.
% 3) For k = 1..N_layers a water-cloud file containing only the top k
%    layers is written; k = 0 is a (practically) cloud-free run. LWC and r_e of every
%    retained layer are IDENTICAL in all files, so R_k - R_(k-1) is exactly
%    the contribution of adding layer k at the base of the cloud.
% 4) w_k = (R_k - R_(k-1)) / dtau_layer / R_N is the layer average of
%    Platnick's eq. (4) over layer k. Summing w_k*dtau_layer over all
%    layers gives 1 - R_0/R_N (R_0 = cloud-free reflectance).
%
% Caveat that ver3 inherits from ver2: write_wc_file computes the layer
% optical depth with its own Mie calculation (gamma distribution,
% inputs.RT.distribution_var), while libRadtran computes optical
% properties from inputs.RT.wc_parameterization ('mie interpolate', the
% libRadtran water-cloud table). If those size distributions differ, the
% optical depth libRadtran actually sees differs slightly from tau_edges.
% I have not quantified this.
% -------------------------------------------------------------------------


clear variables

tic


%% ========================= USER SETTINGS =============================

% --- Cloud definition (Platnick (2000) Figure 1 example: tau_c = 8, 5 -> 12 um) ---
r_top_um = 12;                 % microns - effective radius at cloud top
r_bot_um = 5;                  % microns - effective radius at cloud base
tau_c = 8;                     % total cloud optical depth at inputs.RT.lambda_forTau
profile_type = 'adiabatic';    % passed to create_droplet_profile2

% --- Vertical resolution of the truncation ---
% Optical thickness of each homogeneous layer. N_layers = tau_c/dtau_layer
% truncated clouds (+1 cloud-free run) are computed for every wavelength.
dtau_layer = 0.05;

% Number of uniform-z layers used ONLY to map optical depth to altitude
% for the physical cloud (step 1 above). Each one costs one Mie run.
n_fine_layers = 400;

% --- Smoothing of the weighting functions ---
% Moving-average window (number of layers) applied to w before
% renormalising. 1 = no smoothing. Both the raw w and the smoothed,
% renormalised f are saved.
N_mov_avg = 1;

% --- Optional subset of wavelengths, for quick local tests ---
% Empty = all HySICS bands. Otherwise the band closest to each wavelength
% (nm) is kept, e.g. [650, 1240, 1650, 2130].
wavelengths_subset_nm = [];



%% ======================= Which computer? ===============================

clear inputs

inputs.which_computer = whatComputer;

if strcmp(inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------
    % ver3 uses its own INP / water-cloud folders so it never deletes or
    % overwrites the ver2 working files.

    % Folder where all .INP and .OUT files will be saved
    inputs.folderpath_inp = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/HySICS/weighting_functions_ver3/';

    % libRadtran data path. All paths must be absolute in the INP files.
    inputs.libRadtran_data_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/';

    % Water cloud file location
    inputs.water_cloud_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/wc/weighting_func_ver3/';

    % Mie working folder (write_wc_file writes and runs Mie input files here)
    inputs.mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    % Where the final .mat file is saved
    inputs.folderpath_2save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/weighting_functions/'];

elseif strcmp(inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer --------
    % INP/OUT, water cloud and Mie working files all live on scratch.
    % The libRadtran data path matches ver2.

    inputs.folderpath_inp = '/scratch/alpine/anbu8374/weighting_functions_ver3/';

    inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';

    inputs.water_cloud_folder_path = '/scratch/alpine/anbu8374/weighting_functions_ver3/wc/';

    inputs.mie_folder_path = '/scratch/alpine/anbu8374/weighting_functions_ver3/mie/';

    inputs.folderpath_2save = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/HySICS/weighting_functions/';

else

    error([newline, 'Folder paths are not defined for computer "', inputs.which_computer, '".', newline])

end


% Create any folder that doesn't exist yet
folders_needed = {inputs.folderpath_inp, inputs.water_cloud_folder_path,...
    inputs.mie_folder_path, inputs.folderpath_2save};
for ff = 1:numel(folders_needed)
    if ~exist(folders_needed{ff}, 'dir')
        mkdir(folders_needed{ff})
    end
end



%%  Delete old ver3 working files (only in the ver3 folders)

delete([inputs.folderpath_inp, '*.INP'])
delete([inputs.folderpath_inp, '*.OUT'])
delete([inputs.water_cloud_folder_path, '*.DAT'])



%% Set up the libRadtran inputs (same as ver2)

% Set the number of free parameters for the droplet profile
inputs.RT.num_re_parameters = 2;

% Set the vertical homogeneous string
inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';

% Would you like to print the libRadtran error file?
print_libRadtran_err = true;

% how similar should the forward model be to the simulated measurements?
% options: (1) 'exact'  (2) 'subset'
simulated_measurements_likeness = 'exact';

% set up the inputs to create an INP file for DISORT!
[inputs, spec_response] = create_uvSpec_DISORT_inputs_for_HySICS(inputs, false, [],...
    simulated_measurements_likeness, print_libRadtran_err);

% ----------------------------------------------------------------
% ******************* Redefine a few settings ********************
% ----------------------------------------------------------------

% set the calculation type
inputs.calc_type = 'weighting_functions';

% Black surface. Platnick's interpretation of w_m as a maximum-penetration
% weighting holds for a black surface (Platnick 2000, section 4.1). It also
% guarantees R(tau) cannot decrease when a layer is added at the base.
inputs.RT.surface_albedo = 0;

% compute monochromatic calculations and reflectance
inputs.RT.monochromatic_calc = true;
inputs.RT.compute_reflectivity_uvSpec = true;

% Cloud definition
inputs.RT.r_top = r_top_um;        % microns
inputs.RT.r_bot = r_bot_um;        % microns
inputs.RT.tau_c = tau_c;
inputs.RT.profile_type = profile_type;
inputs.RT.indVar = 'altitude';

% set the sensor altitude at cloud top
inputs.RT.sensor_altitude = inputs.RT.z_topBottom(1);      % km - sensor altitude at cloud top

% Don't modify the total or above-cloud column water vapor
inputs.RT.modify_total_columnWaterVapor = false;
inputs.RT.modify_aboveCloud_columnWaterVapor = false;

% set the geometry (same as ver2)
inputs.RT.sza = 10;           % degree
inputs.RT.phi0 = 0;           % degree
inputs.RT.vza = 19;           % degree
inputs.RT.vaz = 210;          % degree

% Exact backscatter (sza = vza = 0) produced non-monotonic R(tau) in the
% June 2025 truncation run (see header)
if inputs.RT.sza < 1 && inputs.RT.vza < 1
    warning(['sza and vza are both ~0 (exact backscatter). A June 2025 run in this ',...
        'geometry gave R(tau) that decreased as layers were added. Check the ',...
        'negative-w diagnostic below.'])
end

% *** NEVER rescale the cloud optical depth inside libRadtran. ***
% "wc_modify tau set" would rescale the LWC of every retained layer and
% destroy the truncation (see header).
inputs.RT.modify_wc_opticalDepth = false;

% We write the truncated water-cloud files ourselves (section below), so
% write_wc_file's own layer-removal branch is not used.
inputs.compute_weighting_functions = false;

% Every run has a water-cloud file (the k = 0 file has negligible LWC; see Step 3)
inputs.RT.yesCloud = true;



%% Optionally keep only a subset of the HySICS bands

if ~isempty(wavelengths_subset_nm)

    num_bands_all = size(inputs.RT.wavelengths2run, 1);

    % index (into wavelengths2run) of the band closest to each requested wavelength
    idx_keep = zeros(numel(wavelengths_subset_nm), 1);
    for ii = 1:numel(wavelengths_subset_nm)
        [~, idx_keep(ii)] = min(abs(inputs.RT.wavelengths2run(:,1) - wavelengths_subset_nm(ii)));
    end
    idx_keep = unique(idx_keep);

    inputs.bands2run = inputs.bands2run(idx_keep);
    inputs.RT.wavelengths2run = inputs.RT.wavelengths2run(idx_keep, :);

    % Subset every spec_response field whose rows are the bands
    sr_fields = fieldnames(spec_response);
    for ii = 1:numel(sr_fields)
        if size(spec_response.(sr_fields{ii}), 1) == num_bands_all
            spec_response.(sr_fields{ii}) = spec_response.(sr_fields{ii})(idx_keep, :);
        end
    end

end

num_wl = size(inputs.RT.wavelengths2run, 1);



%% Step 1: the physical cloud on a fine uniform-z grid -> tau(z)
% Build the same cloud ver2 builds, but with n_fine_layers thin layers, and
% use write_wc_file's own LWC and extinction to get the optical depth of each
% thin layer. This gives the mapping between altitude and optical depth.

z_top_km = max(inputs.RT.z_topBottom);       % km - cloud top altitude
z_base_km = min(inputs.RT.z_topBottom);      % km - cloud base altitude
H_km = z_top_km - z_base_km;                 % km - geometric thickness

% create_droplet_profile2 normalises z by the min and max of its input, so it
% must be evaluated on a full base-to-top grid; other heights are interpolated.
z_dense_km = linspace(z_base_km, z_top_km, 20001)';      % km - base -> top
re_dense_um = create_droplet_profile2([r_top_um, r_bot_um], z_dense_km,...
    'altitude', profile_type);                           % microns
re_dense_um = re_dense_um(:);

% Fine layers: the same uniform-z edges write_wc_file builds internally
% (linspace from cloud base to cloud top with n_fine_layers+1 levels)
z_fine_edges_km = linspace(z_base_km, z_top_km, n_fine_layers+1)';   % km - base -> top
z_fine_mid_km = (z_fine_edges_km(1:end-1) + z_fine_edges_km(2:end))/2;
re_fine_um = interp1(z_dense_km, re_dense_um, z_fine_mid_km);         % microns - r_e at layer midpoints

% write_wc_file wants one distribution parameter per layer
dist_var_fine = repmat(inputs.RT.distribution_var(1), n_fine_layers, 1);

% Only the LWC (g/m^3, base -> top) and the bulk extinction per unit LWC
% (km^-1 per g/m^3, base -> top) are used; the file it writes is not.
[~, lwc_fine_gm3, ext_fine] = write_wc_file(re_fine_um, tau_c, inputs.RT.z_topBottom,...
    inputs.RT.lambda_forTau, inputs.RT.distribution_str, dist_var_fine,...
    'vert-non-homogeneous', inputs.RT.parameterization_str, 'altitude', false,...
    inputs.which_computer, 1, inputs.RT.num_re_parameters,...
    inputs.water_cloud_folder_path, inputs.mie_folder_path);

% Optical depth of each thin layer: tau = LWC * (extinction per LWC) * dz.
% write_wc_file chooses its LWC slope so these sum to tau_c.
dz_fine_km = H_km / n_fine_layers;
tau_fine_layer = lwc_fine_gm3(:) .* ext_fine(:) .* dz_fine_km;        % base -> top

disp([newline, 'Fine grid: sum of layer optical depths = ', num2str(sum(tau_fine_layer), 6),...
    ' (tau_c = ', num2str(tau_c), ')', newline])

% Cumulative optical depth measured from cloud top (tau = 0) downward,
% and the matching altitudes. LWC is constant within each thin layer, so tau
% is linear in z inside a layer and linear interpolation between edges is exact.
tau_fine_edges_td = [0; cumsum(flipud(tau_fine_layer))];   % top -> base
tau_fine_edges_td = tau_fine_edges_td * tau_c / tau_fine_edges_td(end);   % remove round-off
z_fine_edges_td_km = flipud(z_fine_edges_km);              % top -> base



%% Step 2: re-layer the cloud into N_layers layers of equal optical thickness

N_layers = round(tau_c / dtau_layer);
dtau_layer = tau_c / N_layers;                             % exact layer optical thickness

% Layer edges in optical depth (cloud top = 0) and in altitude, top -> base
tau_edges = linspace(0, tau_c, N_layers+1)';
z_edges_km = interp1(tau_fine_edges_td, z_fine_edges_td_km, tau_edges);
z_edges_km([1, end]) = [z_top_km, z_base_km];

dz_layer_km = -diff(z_edges_km);                           % km - geometric thickness of each layer

% Optical depth at the middle of each layer, and the altitude there
tau_midPoint = (tau_edges(1:end-1) + tau_edges(2:end))/2;  % top -> base
z_mid_km = interp1(tau_fine_edges_td, z_fine_edges_td_km, tau_midPoint);

% Each homogeneous layer gets the effective radius found at its tau-midpoint
% (Platnick (2000, section 5) assigned each layer the size at "the depth
% corresponding to the layer midpoint").
re_layer_um = interp1(z_dense_km, re_dense_um, z_mid_km);  % microns, top -> base

% Exact extinction per unit LWC for these radii, from the same Mie path
% write_wc_file uses. A second call is simpler than reproducing its Mie
% set-up here; again only the extinction output is used.
dist_var_layer = repmat(inputs.RT.distribution_var(1), N_layers, 1);
[~, ~, ext_layer_bu] = write_wc_file(flipud(re_layer_um), tau_c, inputs.RT.z_topBottom,...
    inputs.RT.lambda_forTau, inputs.RT.distribution_str, dist_var_layer,...
    'vert-non-homogeneous', inputs.RT.parameterization_str, 'altitude', false,...
    inputs.which_computer, 2, inputs.RT.num_re_parameters,...
    inputs.water_cloud_folder_path, inputs.mie_folder_path);
ext_layer = flipud(ext_layer_bu(:));                       % km^-1/(g/m^3), top -> base

% LWC that gives every layer exactly dtau_layer: tau = LWC * ext * dz
lwc_layer_gm3 = dtau_layer ./ (ext_layer .* dz_layer_km);  % g/m^3, top -> base

% --- Sanity checks ---
% (a) layers must stay distinct after the 1 mm precision used in the files
if min(dz_layer_km) < 1e-5
    error([newline, 'A layer is thinner than 1 cm; increase dtau_layer.', newline])
end
% (b) layer LWC should be close to the continuous cloud's LWC (linear in
% height from zero at the base) at the layer midpoint. The bottom layer
% spans the steepest relative change in LWC, so it is reported separately.
lwc_slope_gm3_per_km = lwc_fine_gm3(end) / H_km;
lwc_cont_gm3 = lwc_slope_gm3_per_km * (z_mid_km - z_base_km);
rel_diff = lwc_layer_gm3 ./ lwc_cont_gm3 - 1;
disp(['Layer LWC vs continuous LWC at the layer midpoint: max |rel. diff| = ',...
    num2str(max(abs(rel_diff(1:end-1))), 3), ' (all but bottom layer), ',...
    num2str(rel_diff(end), 3), ' (bottom layer)', newline])
disp(['N_layers = ', num2str(N_layers), ', dtau = ', num2str(dtau_layer),...
    ', top-layer dz = ', num2str(1e3*dz_layer_km(1), 3), ' m, bottom-layer dz = ',...
    num2str(1e3*dz_layer_km(end), 3), ' m', newline])

% Keep the inputs structure consistent with the cloud actually simulated
inputs.RT.n_layers = N_layers;
inputs.RT.z_edges = flipud(z_edges_km);                    % km, base -> top (ver2 convention)
inputs.RT.distribution_var = dist_var_layer;



%% Step 3: write one water-cloud file per truncation
% wc_filename{k+1} holds the top k layers.

wc_filename = cell(N_layers+1, 1);

% k = 0 is the cloud-free reference, R_0. It cannot be run with
% inputs.RT.yesCloud = false: write_INP_file writes the wavelength,
% geometry and output lines INSIDE its "if inputs.RT.yesCloud==true" block,
% so a cloud-free INP file would be incomplete. Instead the k = 0 file holds
% the top layer with a negligible LWC (1e-10 g/m^3, optical depth ~1e-11),
% which is cloud-free for all practical purposes.
lwc_negligible_gm3 = 1e-10;
wc_filename{1} = sprintf('WC_ver3_rTop%g_rBot%g_tauC%g_top000_of_%03d_layers.DAT',...
    r_top_um, r_bot_um, tau_c, N_layers);
write_truncated_wc_file(inputs.water_cloud_folder_path, wc_filename{1},...
    z_edges_km(1:2), lwc_negligible_gm3, re_layer_um(1));

for k = 1:N_layers

    wc_filename{k+1} = sprintf('WC_ver3_rTop%g_rBot%g_tauC%g_top%03d_of_%03d_layers.DAT',...
        r_top_um, r_bot_um, tau_c, k, N_layers);

    write_truncated_wc_file(inputs.water_cloud_folder_path, wc_filename{k+1},...
        z_edges_km(1:k+1), lwc_layer_gm3(1:k), re_layer_um(1:k));

end



%% Step 4: write the INP files and run DISORT
% One run per (wavelength, truncation). k varies fastest.

num_states = N_layers + 1;                     % k = 0 ... N_layers
num_INP_files = num_states * num_wl;

k_all = repmat((0:N_layers)', num_wl, 1);                          % layers kept
wl_idx_all = reshape(repmat(1:num_wl, num_states, 1), [], 1);      % row of wavelengths2run

% Variables by column:
% (1) - monochromatic wavelength to run (nm)
% (2) - optical depth of the base of the truncated cloud, tau_k
% (3) - k, the number of layers kept (0 = cloud-free)
% (4) - index of the band in inputs.RT.wavelengths2run / spec_response
changing_variables = [inputs.RT.wavelengths2run(wl_idx_all, 1), tau_edges(k_all+1),...
    k_all, wl_idx_all];

inputFileName = cell(num_INP_files, 1);
outputFileName = cell(num_INP_files, 1);
for nn = 1:num_INP_files
    inputFileName{nn} = sprintf('wf_ver3_%.3fnm_k%03d.INP', changing_variables(nn,1), k_all(nn));
    outputFileName{nn} = ['OUTPUT_', inputFileName{nn}(1:end-4)];
end

% store the reflectances
Refl_model = zeros(num_INP_files, 1);

parfor nn = 1:num_INP_files

    % ------------------ Write the INP File --------------------
    % No wc_modify argument: the cloud optical depth is never rescaled.
    % k = 0 uses the negligible-LWC file (see Step 3).
    write_INP_file(inputs.folderpath_inp, inputs.libRadtran_data_path, inputs.water_cloud_folder_path,...
        inputFileName{nn}, inputs, changing_variables(nn,1), wc_filename{k_all(nn)+1});

    % ----------------- Run radiative transfer -----------------
    runUVSPEC_ver2(inputs.folderpath_inp, inputFileName{nn}, outputFileName{nn},...
        inputs.which_computer);

    % read .OUT file
    [ds,~,~] = readUVSPEC_ver2(inputs.folderpath_inp, outputFileName{nn}, inputs,...
        inputs.RT.compute_reflectivity_uvSpec);

    % Store the azimuthally averaged reflectance (1/sr)
    Refl_model(nn) = ds.reflectivity.az_avg;

end

toc



%% Step 5: weighting functions, Platnick (2000) eq. 4

% Rows: k = 0 (cloud-free), 1, ..., N_layers. Columns: wavelengths.
Refl_model = reshape(Refl_model, num_states, num_wl);

R_clear = Refl_model(1, :);          % cloud-free reflectance
R_full = Refl_model(end, :);         % reflectance of the full tau_c cloud

% Layer-average of dR/dtau for the layer added at the base, normalised by
% R(tau_c). Size: [N_layers x num_wl], top -> base.
w = diff(Refl_model, 1, 1) ./ dtau_layer ./ R_full;

% --- Diagnostics ---
% Summed over layers, w integrates to 1 - R_clear/R_full (eq. 4 telescopes)
area_w = sum(w, 1) * dtau_layer;
disp(['Max |sum(w*dtau) - (1 - R_clear/R_full)| = ',...
    num2str(max(abs(area_w - (1 - R_clear./R_full))), 3), newline])

% Over a black surface R(tau) cannot decrease when a layer is added at the
% base, so any negative w is numerical noise. Report how large it is
% relative to the peak of that band's weighting function.
num_neg = nnz(w < 0);
rel_neg = min(w, [], 1) ./ max(w, [], 1);     % per band
[worst_rel_neg, idx_worst] = min(rel_neg);
disp([num2str(num_neg), ' of ', num2str(numel(w)), ' weighting-function values are negative.'])
if num_neg > 0
    disp(['Most negative value relative to its band peak: ', num2str(worst_rel_neg, 3),...
        ' at ', num2str(inputs.RT.wavelengths2run(idx_worst,1)), ' nm', newline])
end



%% Optional smoothing, then renormalise to unit area
% f is what the information-content analysis uses. Renormalising makes eq. 3
% return r_e* = r_e for a vertically uniform profile (Platnick 2000, section 4).

f = w;
if N_mov_avg > 1
    f = movmean(w, N_mov_avg, 1);
end
f = f ./ (sum(f, 1) * dtau_layer);           % layer quadrature: sum(f*dtau) = 1



%% Save the output

rev = 1;
filename_base = [inputs.folderpath_2save, 'disort_HySICS_reflectance_for_weightingFunctions_',...
    'truncatedCloud_ver3_', profile_type, '_rTop', num2str(r_top_um), '_rBot', num2str(r_bot_um),...
    '_tauC', num2str(tau_c), '_sim-ran-on-', char(datetime("today"))];
filename = [filename_base, '_rev', num2str(rev), '.mat'];
while isfile(filename)
    rev = rev + 1;
    filename = [filename_base, '_rev', num2str(rev), '.mat'];
end

weighting_function_method = ['Platnick (2000) eq. 4: R(tau) from truncated clouds ',...
    '(lower layers absent), equal-tau homogeneous layers, black surface'];

save(filename, "Refl_model", "R_clear", "R_full", "inputs", "spec_response",...
    "changing_variables", "w", "f", "N_mov_avg", "tau_edges", "tau_midPoint",...
    "dtau_layer", "N_layers", "re_layer_um", "lwc_layer_gm3", "z_edges_km",...
    "weighting_function_method");

disp(['Saved: ', filename, newline])



%% Plot a few weighting functions and R(tau)

wl_2plot_nm = [650, 1240, 1650, 2130];       % nearest computed bands are used
idx2plot = zeros(numel(wl_2plot_nm), 1);
for ii = 1:numel(wl_2plot_nm)
    [~, idx2plot(ii)] = min(abs(inputs.RT.wavelengths2run(:,1) - wl_2plot_nm(ii)));
end
idx2plot = unique(idx2plot);

lgnd_str = cell(numel(idx2plot), 1);
for ii = 1:numel(idx2plot)
    lgnd_str{ii} = ['$\lambda = $', num2str(round(inputs.RT.wavelengths2run(idx2plot(ii), 1), 1)), ' $nm$'];
end

figure;

% Weighting functions (raw w; f is identical when N_mov_avg = 1)
subplot(1,2,1)
for ii = 1:numel(idx2plot)
    plot(w(:, idx2plot(ii)), tau_midPoint, 'Color', mySavedColors(ii, 'fixed'))
    hold on
end
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$w_m(\tau)$', 'Interpreter', 'latex')
ylabel('$\tau$', 'Interpreter', 'latex')
title('Truncated-cloud weighting functions (ver3)', 'Interpreter', 'latex')
legend(lgnd_str, 'Interpreter', 'latex', 'Location', 'southeast')

% Reflectance of the truncated cloud versus the optical depth of its base
subplot(1,2,2)
for ii = 1:numel(idx2plot)
    plot(Refl_model(:, idx2plot(ii)), tau_edges, 'Color', mySavedColors(ii, 'fixed'))
    hold on
end
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$R(\tau)$ (1/sr)', 'Interpreter', 'latex')
ylabel('optical depth of the truncated cloud base', 'Interpreter', 'latex')
title(['$r_{top} = $', num2str(r_top_um), ' $\mu m$, $r_{bot} = $', num2str(r_bot_um),...
    ' $\mu m$, $\tau_c = $', num2str(tau_c)], 'Interpreter', 'latex')

set(gcf, 'Position', [0 0 1400 800])



%% ======================== Local functions ==============================

function write_truncated_wc_file(folder_path, file_name, z_edges_td_km, lwc_td_gm3, re_td_um)
% Write a libRadtran "wc_file 1D" holding only the top layers of the cloud.
%
% INPUTS (all ordered from cloud top downward):
%   z_edges_td_km - layer edges (km above sea level), length = number of layers + 1
%   lwc_td_gm3    - liquid water content of each layer (g/m^3)
%   re_td_um      - effective radius of each layer (microns)
%
% libRadtran (>= 1.4) interprets the values on a line as properties of the
% LAYER ABOVE that altitude (libRadtran 2.0.6 doc/radiative_transfer.tex,
% subsection "Water clouds": with lines at 5, 4, 3, 2 km "the cloud would
% extend from 2 to 5 km, with e.g. a LWC of 0.2 g/m^3 for the layer between
% 4 and 5 km"). write_wc_file uses the same convention. So:
%   - the top line (cloud-top altitude) carries zeros: no cloud above it;
%   - the line at the base of layer i carries layer i's LWC and r_e;
%   - a final line at 0 km with zeros makes the air between the surface and
%     the truncated cloud base cloud-free.
% Altitudes are written with 1 mm precision and LWC with 1e-10 g/m^3
% precision. The %12.3f / %7.4f format in write_wc_file would merge the
% metre-thin layers near cloud top and round the small LWCs near the base.

fileID = fopen([folder_path, file_name], 'w');
if fileID < 0
    error([newline, 'Could not open ', folder_path, file_name, ' for writing.', newline])
end

fprintf(fileID, '%s %13s %16s %12s\n', '#', 'z', 'LWC', 'R_eff');
fprintf(fileID, '%s %13s %16s %12s\n', '#', '(km)', '(g/m^3)', '(micron)');

% cloud top: nothing above
fprintf(fileID, '%15.6f %16.10f %12.6f\n', z_edges_td_km(1), 0, 0);

% each retained layer, written at its lower edge
for ii = 1:numel(lwc_td_gm3)
    fprintf(fileID, '%15.6f %16.10f %12.6f\n', z_edges_td_km(ii+1), lwc_td_gm3(ii), re_td_um(ii));
end

% clear air between the surface and the truncated cloud base
if z_edges_td_km(end) > 0
    fprintf(fileID, '%15.6f %16.10f %12.6f\n', 0, 0, 0);
end

fclose(fileID);

end
