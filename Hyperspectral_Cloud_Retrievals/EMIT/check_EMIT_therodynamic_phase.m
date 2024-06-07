%% Check the thermodynamic phase of some EMIT pixel

% Thermodynamic phase discrimination is based off of the method outlined by
% D.R. Thompson et al. (2016): "Measuring cloud thermodynamic phase with
% shortwave infrared imaging spectroscopy"

% ------------------------------------------------------------------
% ----------------------- THINGS TO WORK ON ------------------------
% ------------------------------------------------------------------
% This function is not yet properly working. It can compute the bulk
% absorption coefficient for all three phases of water, and it can solve
% the non-negative least squares problem. But what how do I take this
% information and determine the phase of the cloud? What if it's mixed
% phase? I haven't verified the results yet. I can do this by running
% Yolanda's code to find an EMIT data set close to some other satellite
% instrument like GOES, VIIRS, or MODIS, and verify the phase I compute
% with the ansewr these instruments find. Do I need to add the bulk
% absorption coefficient of other absorbers? From figure 2.12 in Bohren and
% Clothiaux's book, Fundamentals of Atmospheric Radiation, there are a few
% other species that I might want to model: CO2, CO, Methane and Nitrous
% Oxide.

% --------------------------------------------------------------------
% --------------------------------------------------------------------


% By Andrew John Buggee

function inputs = check_EMIT_therodynamic_phase(emit, inputs)


%% We start by assuming if we deal with small wavelength intervals, the reflectance is linear with lambda
% We will use the region from 1400 to 1800 nms

wavelength_boundaries = [1400, 1800];       % nm
%wavelength_boundaries = [300, 2500];       % nm


% To make our linear assumption valid, let's use the wavelength spacing
% between the center wavelengths of adjacent emit channels
wavelength_idx = emit.radiance.wavelength >= wavelength_boundaries(1)...
    & emit.radiance.wavelength <= wavelength_boundaries(2);

reflectance = emit.reflectance.value(wavelength_idx, :);      % 1/sr

% the center wavelengths for each measurement
wavelength_center_emit = emit.radiance.wavelength(wavelength_idx);       % nm


% For our calculations of bulk absorption coefficients, let's use a high
% resolution wavelength grid that can be down sampled later to match the
% wavelength grid of EMIT measurements. Make sure the grid is bounded by
% the wavelengths on the EMIT wavelength grid
grid_resolution = 1;        % nm
wavelength_grid_fine = wavelength_center_emit(1):grid_resolution:wavelength_center_emit(end);  % nm

%% Using Mie theory, compute the absorption coefficients

% We must do this for all atmospheric absorbers
% ------------------------------------------------------------------
% Run Mie calculation to obtain single scatering albedo and asymmetry
% parameter
% ------------------------------------------------------------------

% The first entry belows describes the type of droplet distribution
% that should be used. The second describes the distribution width. If
% running a mono-dispersed calculation, no entry for distribution width is
% required.
mie_distribution = {'mono', []};           % droplet distribution

% What mie code should we use to compute the scattering properties?
mie_program = 'MIEV0';               % type of mie algorithm to run

% Do you want a long or short error file?
err_msg_str = 'verbose';

% What is the index of refraction of the scatterer? If there is more than
% one scatterer, enter multiple values for the index of refraction
indexOfRefraction = {'water', 'ice'};

% Define the size of the scatterer and its scattering properties
% Assuming a pure homogenous medium composed of a single substance.
% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
radius = 10;            % microns

% **** DON'T NEED FINE GRID FOR LIBRADTRAN COMPUTATIONS ****
% we don't need to use a finer grid for liquid and solid phases of
% water. There are now downsampling issues.

% Step through each absorber and each wavelength
k_bulk_water = zeros(1, length(wavelength_center_emit));
k_bulk_ice = zeros(1, length(wavelength_center_emit));

for jj = 1:length(indexOfRefraction)


    for ww = 1:length(wavelength_center_emit)

        disp(['jj=',num2str(jj),', ww=',num2str(ww)])
        % define the wavelength
        % The wavelength input is defined as follows:
        % [wavelength_start, wavelength_end, wavelength_step].
        wl = [wavelength_center_emit(ww), wavelength_center_emit(ww), 0];

        % Create a mie file
        [input_filename, output_filename, mie_folder] = write_mie_file(mie_program, indexOfRefraction{jj},...
            radius, wl, mie_distribution, err_msg_str, 1);

        % run the mie file
        runMIE(mie_folder,input_filename,output_filename);

        % Read the output of the mie file
        [ds,~,~] = readMIE(mie_folder,output_filename);

        % store the complex index of refraction
        if jj==1

            k_bulk_water(ww) = 4*pi*ds.refrac_imag/(wavelength_center_emit(ww)*1e-7);       % centimeters^(-1)

        elseif jj==2

            k_bulk_ice(ww) = 4*pi*ds.refrac_imag/(wavelength_center_emit(ww)*1e-7);       % centimeters^(-1)

        end

    end

end

%% Compute the bulk absorption coefficient for water vapor

hitran_waterVapor_file = 'hitran_water_vapor_absorption_350_to_2600nm.mat';

% --- Define the pressure, temperature and wavelength range ---
con = physical_constants();
% read in standard values of atmospheric pressure and temperature at some
% reasonable height for the lifting condensation level using the US
% standard atmosphere

H_lcl = 3;            % km
atm_prof = read_atm_profile(H_lcl, 'atm_profile_with_homogenous_cloud.txt', true);

T = atm_prof.temperature;                % K - temperature of gas
P = (atm_prof.pressure*100)/con.atm;                  % atm - pressure of whole atmosphere

% because number density is directly proportional to the pressure, the
% partial pressure is the ratio of number density of water vapor to the
% number density of air.
P_self = atm_prof.H2O_Nc/atm_prof.air_Nc;          % atm - partial pressure of water vapor


% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

tic
abs_waterVapor = hitran_compute_abs_cross_section(hitran_waterVapor_file, T, P, P_self,...
    wavelength_grid_fine, solution_type);

disp([newline,'Time it took to compute the bulk absorption coefficient for ',...
    'water vapor: ', num2str(toc), ' sec', newline])

%% Compute the bulk absorption coefficient for carbon dioxide

hitran_carbonDioxide_file = 'hitran_carbonDioxide_absorption_350_to_2600nm.mat';

% --- Define the pressure, temperature and wavelength range ---
con = physical_constants();
% read in standard values of atmospheric pressure and temperature at some
% reasonable height for the lifting condensation level using the US
% standard atmosphere

H_lcl = 3;            % km
atm_prof = read_atm_profile(H_lcl, 'atm_profile_with_homogenous_cloud.txt', true);

T = atm_prof.temperature;                % K - temperature of gas
P = (atm_prof.pressure*100)/con.atm;                  % atm - pressure of whole atmosphere

% because number density is directly proportional to the pressure, the
% partial pressure is the ratio of number density of Carbon Dioxide to the
% number density of air.
P_self = atm_prof.CO2_Nc/atm_prof.air_Nc;          % atm - partial pressure of water vapor

% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

abs_co2 = hitran_compute_abs_cross_section(hitran_carbonDioxide_file, T, P, P_self,...
    wavelength_grid_fine, solution_type);


%% Compute the bulk absorption coefficient of Methane

hitran_methane_file = 'hitran_methane_absorption_350_to_2600nm.mat';

% --- Define the pressure, temperature and wavelength range ---
con = physical_constants();
% read in standard values of atmospheric pressure and temperature at some
% reasonable height for the lifting condensation level using the US
% standard atmosphere

H_lcl = 3;            % km
atm_prof = read_atm_profile(H_lcl, 'atm_profile_with_homogenous_cloud.txt', true);

T = atm_prof.temperature;                % K - temperature of gas
P = (atm_prof.pressure*100)/con.atm;                  % atm - pressure of whole atmosphere

% Partial pressure of methane from values listed at:
% https://www.e-education.psu.edu/meteo300
P_self = 0.00000182 * P;          % atm - partial pressure of methane

% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

abs_methane = hitran_compute_abs_cross_section(hitran_methane_file, T, P, P_self,...
    wavelength_grid_fine, solution_type);


%% Compute the bulk absorption coefficient of Carbon Monoxide

hitran_carbonMonoxide_file = 'hitran_carbonMonoxide_absorption_350_to_2600nm.mat';

% --- Define the pressure, temperature and wavelength range ---
con = physical_constants();
% read in standard values of atmospheric pressure and temperature at some
% reasonable height for the lifting condensation level using the US
% standard atmosphere

H_lcl = 3;            % km
atm_prof = read_atm_profile(H_lcl, 'atm_profile_with_homogenous_cloud.txt', true);

T = atm_prof.temperature;                % K - temperature of gas
P = (atm_prof.pressure*100)/con.atm;                  % atm - pressure of whole atmosphere

% Partial pressure of carbon monoxide from values listed at:
% https://scied.ucar.edu/learning-zone/air-quality/carbon-monoxide
P_self = 100e-9 * P;          % atm - partial pressure of methane

% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

abs_CO = hitran_compute_abs_cross_section(hitran_carbonMonoxide_file, T, P, P_self,...
    wavelength_grid_fine, solution_type);


%% Compute the bulk absorption coefficient of molecular oxygen

hitran_oxygen_file = 'hitran_oxygen_absorption_350_to_2600nm.mat';

% --- Define the pressure, temperature and wavelength range ---
con = physical_constants();
% read in standard values of atmospheric pressure and temperature at some
% reasonable height for the lifting condensation level using the US
% standard atmosphere

H_lcl = 3;            % km
atm_prof = read_atm_profile(H_lcl, 'atm_profile_with_homogenous_cloud.txt', true);

T = atm_prof.temperature;                % K - temperature of gas
P = (atm_prof.pressure*100)/con.atm;                  % atm - pressure of whole atmosphere

% Partial pressure of carbon monoxide from values listed at:
% https://scied.ucar.edu/learning-zone/air-quality/carbon-monoxide
P_self = 100e-9 * P;          % atm - partial pressure of methane

% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

abs_O2 = hitran_compute_abs_cross_section(hitran_oxygen_file, T, P, P_self,...
    wavelength_grid_fine, solution_type);



%% Set up the non-negative least squares solution

% Following D.R. Thompson et al. 2016

% Solve the problem ||d - Gm|| subject to m>0



% set up the matrix G

% *** QUESTION ***
% Should I take the average value of the bulk absorption coeff over each
% EMIT channel?

% G = [ones(length(wavelength_center_emit), 1), wavelength_center_emit, -wavelength_center_emit,...
%     k_bulk_water', k_bulk_ice', abs_waterVapor.bulk_coefficient, abs_co2.bulk_coefficient,...
%     abs_methane.bulk_coefficient, abs_CO.bulk_coefficient, abs_O2.bulk_coefficient];
%
% % Compute the non-negative least squares solution for each pixel
% num_pix = size(reflectance,2);
% x = zeros(size(G,2), num_pix);
%
% % Make sure to translate the reflectance measurements into negative log
% % space!
%
% for pp = 1:num_pix
%
%     x(:,pp) = lsqnonneg(G, -log(reflectance(:, pp)));
%
% end


% the bulk absorption coefficients for liquid water and ice are smooth
% enough that we can linear interpolate to get it in on the same fine grid.
k_bulk_water_fine = interp1(wavelength_center_emit, k_bulk_water, wavelength_grid_fine);
k_bulk_ice_fine = interp1(wavelength_center_emit, k_bulk_ice, wavelength_grid_fine);

% Create the G matrix
G_fine = [ones(length(wavelength_grid_fine), 1), wavelength_grid_fine', -wavelength_grid_fine',...
    k_bulk_water_fine', k_bulk_ice_fine', abs_waterVapor.bulk_coefficient, abs_co2.bulk_coefficient,...
    abs_methane.bulk_coefficient, abs_CO.bulk_coefficient, abs_O2.bulk_coefficient];

% Compute the non-negative least squares solution for each pixel
num_pix = size(reflectance,2);
x = zeros(size(G_fine,2), num_pix);
reflectance_fine = zeros(length(wavelength_grid_fine), num_pix);


% Make sure to translate the reflectance measurements into negative log
% space!

for pp = 1:num_pix

    % Map the reflectance measurements onto the finer grid. There are some sharp
    % deviations, but only a few and linear interpolation is able to capture
    % them
    reflectance_fine(:,pp) = interp1(wavelength_center_emit, reflectance(:,pp), wavelength_grid_fine);

    % compute the non-negative least squares solution
    x(:,pp) = lsqnonneg(G_fine, -log(reflectance_fine(:, pp)));

end


% Compute the simulated reflectances using the solution to the non-negative
% least squares problem
reflectance_sim = G_fine*x;
reflectance_sim = exp(-reflectance_sim);

% --- Store the retrieved values of the effective water thickness ---
inputs.phase.effectiveWaterThickness.liquid = x(4,:);      % liquid water thickness in cm
inputs.phase.effectiveWaterThickness.ice = x(5,:);      % ice water thickness in cm
inputs.phase.effectiveWaterThickness.vapor = x(6,:);      % water vapor thickness in cm

% also store the simulated reflectance and the wavelength grid
inputs.phase.reflectance_sim = reflectance_sim;
inputs.phase.wavelength_grid = wavelength_grid_fine;

%% Compare the true reflectance with the simulated reflectance

figure;
for nn = 1:num_pix
    plot(wavelength_center_emit, reflectance(:,nn),'.-', 'markersize', 20,...
        'linewidth',1, 'Color', mySavedColors(nn,'fixed'))
    hold on;
    plot(wavelength_grid_fine, reflectance_sim(:,nn), '-', 'Color',...
        mySavedColors(nn,'fixed'))


end

grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('Reflectance (1/sr)', 'Interpreter','latex')
legend('Pixel 1 - EMIT', 'Pixel 1 - Calculated', 'Pixel 2 - EMIT', 'Pixel 2 - Calculated', ...
    'Location', 'best')
set(gcf, 'Position', [0 0 900 400])

title('$H_{2}O$ / $CO_{2}$ / $CH_{4}$ / $O_{2}$ / $CO$', 'Interpreter','latex')

end

