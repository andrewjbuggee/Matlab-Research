%% Check the thermodynamic phase of some EMIT pixel

% Thermodynamic phase discrimination is based off of the method outlined by
% D.R. Thompson et al. (2016): "Measuring cloud thermodynamic phase with
% shortwave infrared imaging spectroscopy"


% By Andrew John Buggee

function inputs = check_EMIT_therodynamic_phase(emit, pixels2use)


%% We start by assuming if we deal with small wavelength intervals, the reflectance is linear with lambda
% We will use the region from 1400 to 1800 nms

wavelength_boundaries = [1400, 1800];       % nm

% To make out linear assumption valid, let's use the wavelength spacing
% between the center wavelengths of adjacent emit channels
wavelength_idx = emit.radiance.wavelength >= wavelength_boundaries(1)...
    & emit.radiance.wavelength <= wavelength_boundaries(2);

reflectance = emit.reflectance(wavelength_idx, :);      % 1/sr

% the center wavelengths for each measurement
%wavelength_center_emit = emit.radiance.wavelength(wavelength_idx);       % nm
wavelength_center_emit = (wavelength_boundaries(1):5:wavelength_boundaries(2))';        % nm


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


% Step through each absorber and each wavelength
k_bulk_water = zeros(1, length(wavelength_center_emit));
k_bulk_ice = zeros(1, length(wavelength_center_emit));

for jj = 1:length(indexOfRefraction)


    for ww = 1:length(wavelength_center_emit)

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

% --- Define the pressure, temperature and wavelength range
T = 298;                % K - temperature of gas
P = 1;                  % atm - pressure of whole atmosphere
P_self = 0.03;          % atm - partial pressure of water vapor

% --- Define the solution type ----
% How should we solve for the absorption cross section?
% 'full_integral' solves for the Voigt function by solving a numerical
% integral over a finely spaced grid. This takes a while
% 'whitting' uses an approximation
solution_type = 'whitting';

abs_waterVapor = hitran_compute_abs_cross_section(hitran_waterVapor_file, T, P, P_self,...
    wavelength_center_emit, solution_type);


%% Set up the non-negative least squares solution

% Following D.R. Thompson et al. 2016

% Solve the problem ||d - Gm|| subject to m>0

% set up the matrix G

% *** QUESTION ***
% Should I take the average value of the bulk absorption coeff over each
% EMIT channel?

G = [ones(length(wavelength_center_emit), 1), wavelength_center_emit', -wavelength_center_emit',...
    k_bulk_water', k_bulk_ice', abs_waterVapor.bulk_coefficient];





end

