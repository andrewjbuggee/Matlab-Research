% simulated measurements is an input used to define certain parameters

% ** Retrieving 4 variables: log(r_top), log(r_bot), log(tau_c), log(cwv)


function GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS_ver4_logState(simulated_measurements, print_libRadtran_err)


% Which computer are you using?
GN_inputs.which_computer = whatComputer;





%%

% define the number of iterations for the gauss-newton solver
GN_inputs.GN_iterations = 5;


% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
GN_inputs.percent_change_limit = 0.03;

% define the type of model prior pdf
GN_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
GN_inputs.num_model_parameters = 4;




% -------------------------------------------
% --- Stuff for the Model Parameter Prior ---
% -------------------------------------------

    

GN_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth', 'Above Cloud Column Water Vapor'};


% ---------------------------------------
% --- Stuff for the Measurement Prior ---
% ---------------------------------------


GN_inputs.measurement.prior = 'gaussian';
% covaraince_type can be:
%   (1) 'independent - thus all off diagonal elements are 0
%   (2) 'computed' - uses measured data to compute covaraince
GN_inputs.measurement.covariance_type = 'independent';

% Store the simulated state vector used to create the measurements
if isfield(simulated_measurements.inputs.RT, 'r_top')
    % Then we have a simulated adiabatic or other theoretical profile
    GN_inputs.measurement.r_top = simulated_measurements.inputs.RT.r_top;      % microns
    GN_inputs.measurement.r_bot = simulated_measurements.inputs.RT.r_bot;      % microns

elseif isfield(simulated_measurements.inputs.RT, 're')

    GN_inputs.measurement.re_prof = simulated_measurements.inputs.RT.re;   % in-situ re profile

    if isfield(simulated_measurements.inputs.RT, 'tau')==true

        GN_inputs.measurement.tau_prof = simulated_measurements.inputs.RT.tau;   % in-situ dervied optical depth vector

    end

    GN_inputs.measurement.lwc_prof = simulated_measurements.inputs.RT.lwc;   % in-situ lwc profile
    GN_inputs.measurement.z = simulated_measurements.inputs.RT.z;            % altidue vector for in-situ measurements
end

GN_inputs.measurement.tau_c = simulated_measurements.inputs.RT.tau_c;      % optical depth
GN_inputs.measurement.actpw = aboveCloud_CWV_simulated_hysics_spectra(simulated_measurements.inputs); % kg/m^2 (equivelant to mm)


% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

GN_inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';

% we model two free parameters, r_top and r_bot
GN_inputs.RT.num_re_parameters = 2;

% we have to assume a vertical profile of droplet size with cloud optical
% depth exsists. And we only retrieve the droplet size at the top and
% bottom. This is the method of King and Vaughn (2012)

% the options for the vertical droplet profile are:
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
% x is determined by the choice of droplet profile within the function
% create_droplet_profile.m

GN_inputs.model.profile.type = 'adiabatic';
GN_inputs.model.profile.r_top = 10; % microns - value for our model
GN_inputs.model.profile.r_bottom = 5; % microns - value for our model







% -----------------------------------------------
% --------------- Define flags  -----------------
% -----------------------------------------------










% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------

% Define the parameters of the INP file

% Use geometry inputs from the simulated measurements
load_parameters_from_measurement = true;

% how similar should the forward model be to the simulated measurements?
% options: (1) 'exact'  (2) 'subset'

simulated_measurements_likeness = 'exact';

[GN_inputs, ~] = create_uvSpec_DISORT_inputs_for_HySICS(GN_inputs, load_parameters_from_measurement,...
    simulated_measurements, simulated_measurements_likeness, print_libRadtran_err);

% Are you simulating a measurement, or making forward model calculations
% for the retrieval?
GN_inputs.calc_type = 'forward_model_calcs_forRetrieval';











end
