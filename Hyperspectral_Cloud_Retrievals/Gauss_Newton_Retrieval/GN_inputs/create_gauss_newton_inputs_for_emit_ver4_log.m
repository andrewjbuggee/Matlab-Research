%% ----- CREATE INPUTS NEEDED TO COMPUTE DROPLET RADIUS PROFILE USING HYPERSPECTRAL MEASUREMENTS FROM EMIT -----

% ** Retrieving 4 variables: log(r_top), log(r_bot), log(tau_c), log(cwv)

% INPUTS:
%   (1) emitDataFolder - 

%   (2) folder2save - 

%   (3) L1B_filename - 

%   (4) emit - 


% OUTPUTS:
%   (1) GN_inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function GN_inputs = create_gauss_newton_inputs_for_emit_ver4_log(emit, print_libRadtran_err)


% Which computer are you using?
GN_inputs.which_computer = whatComputer;



%%

 

% -----------------------------------------------------------
% --------------- Iterations and Convergence ----------------
% -----------------------------------------------------------

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
% --------------- Define flags  -----------------
% -----------------------------------------------
% define flags that tell the codes to either run certain things, or don't
% run certain things

GN_inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use

% if true, the code will load an older set of pixels that has already been used before, and 
% likely has INP files. If false, it tells the code to find a new random subset of pixels
GN_inputs.flags.loadPixelSet = true; 
GN_inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
GN_inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
GN_inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot 





% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------

% Define the parameters of the INP file

% define whether this is a vertically homogenous cloud or not
GN_inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';


GN_inputs = create_uvSpec_DISORT_inputs_for_EMIT(GN_inputs, emit, print_libRadtran_err);










% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check GN_inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
