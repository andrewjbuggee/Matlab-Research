% simulated measurements is an input used to define certain parameters

% ** Retrieving 4 variables: r_top, r_bot, tau_c, cwv


function GN_inputs = create_gauss_newton_inputs_for_simulated_HySICS_ver2(simulated_measurements)


% Which computer are you using?
GN_inputs.which_computer = whatComputer;


% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(GN_inputs.which_computer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    GN_inputs.libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';

   


elseif strcmp(GN_inputs.which_computer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    GN_inputs.libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];

   


elseif strcmp(GN_inputs.which_computer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------


    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    GN_inputs.libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';



end




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

% Using the King and Vaughn (2012) method, we retireve 3 parameters
%   (1) r_top = effective droplet size at the cloud top
%   (2) r_bottom = effective droplet size at the cloud bottom
%   (3) tau_c = cloud optical depth
% a good starting place is to assume the droplet size at cloud top and
% bottom are the same value



    

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

% -----------------------------------------------
% --- Stuff for the Assumed Vertical Profile ---
% -----------------------------------------------

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
% ------------- Folder Locations  ---------------
% -----------------------------------------------
GN_inputs.save_calcs_fileName = ['uvspec_GaussNewton_calcs_',date,'.mat'];




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

[GN_inputs, ~] = create_uvSpec_DISORT_inputs_for_HySICS(GN_inputs, load_parameters_from_measurement, ...
    simulated_measurements, simulated_measurements_likeness);










end
