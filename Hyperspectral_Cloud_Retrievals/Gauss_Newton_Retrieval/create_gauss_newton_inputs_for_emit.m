%% ----- CREATE INPUTS NEEDED TO COMPUTE DROPLET RADIUS PROFILE USING HYPERSPECTRAL MEASUREMENTS FROM EMIT -----


% INPUTS:
%   (1) emitDataFolder - 

%   (2) folder2save - 

%   (3) L1B_filename - 

%   (4) emit - 


% OUTPUTS:
%   (1) GN_inputs - effective droplet radius profile


% By Andrew John Buggee
%%

function GN_inputs = create_gauss_newton_inputs_for_emit(emitDataFolder, folderpaths, L1B_fileName, emit)


% Which computer are you using?
GN_inputs.which_computer = whatComputer;


% --- SAVE THE EMIT FILE NAME ----
GN_inputs.emitDataFolder = emitDataFolder;


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

% ----- Save the L1B file name -----
GN_inputs.L1B_filename = L1B_fileName{1};
 

% -----------------------------------------------------------
% --------------- Iterations and Convergence ----------------
% -----------------------------------------------------------

% define the number of iterations for the gauss-newton solver
GN_inputs.GN_iterations = 5;


% define a percent threshold of the difference between successive
% iterations. If the percent difference is below the percent threshold,
% than the iterative process is stopped.
GN_inputs.percent_change_limit = 0.03;


% Define the convergence limit. Convergence is defined using the residual,
% which is the difference between the true and estimated measurements.
% We take the RMS of the residual using all spectral channels. This is how
% we define the convergence limit. If the residual is the difference
% between the true measurement and the estimated measurement, and the true
% measurement has an uncertainty of 10%, then our estimate measurement
% should be within this uncertainty. Using MODIS, we can compute the
% RMS uncertainty vector and set this as the convergence limit.

GN_inputs.convergence_limit = 0.005;  % generic convergence limit


% define the type of model prior pdf
GN_inputs.model.prior = 'gaussian';


% define the number of model parameters to solve for
GN_inputs.num_model_parameters = 3;


GN_inputs.model.param_names = {'Effective Radius at Top of Cloud', 'Effective Radius at Bottom of Cloud',...
    'Cloud Optical Depth'};



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
% define flags that tell the codes to either run certain things, or don't
% run certain things

GN_inputs.flags.findSuitablePixels = false; % if true, this will search the modis data set for pixels to use

% if true, the code will load an older set of pixels that has already been used before, and 
% likely has INP files. If false, it tells the code to find a new random subset of pixels
GN_inputs.flags.loadPixelSet = true; 
GN_inputs.flags.writeINPfiles = true; % if true, this will create inp files for each the length of vector pixel.row
GN_inputs.flags.runUVSPEC = true; % if true, this will run all of the inp files create from the above flag through uvspec
GN_inputs.flags.plotMLS_figures = false; % this will tell the leasSquaresGridSearch code to plot 





% --------------------------------------------
% Create a new folder to save all calculations
% --------------------------------------------

% Define the folder that stores the GN_inputs and calculated reflectanes
% using todays date
data_date = datetime([L1B_fileName{1}(18:21), '-', L1B_fileName{1}(22:23), '-', L1B_fileName{1}(24:25)],...
    'InputFormat','yyyy-MM-dd');

% Store the file name for the libRadTran INP and OUT files
GN_inputs.folder2save.libRadtran_inp = [folderpaths.libRadtran_inp, 'EMIT_',char(data_date),...
    '_time_', L1B_fileName{1}(27:30), '/'];


% This is the folder where the reflectance calculations will be stored
GN_inputs.folder2save.reflectance_calcs = [folderpaths.reflectance_calcs, emitDataFolder]; 

% If the folder path doesn't exit, create a new directory
if ~exist(GN_inputs.folder2save.reflectance_calcs, 'dir')

    mkdir(GN_inputs.folder2save.reflectance_calcs);

end

% This is the name of the .mat file with the reflectance calcs
GN_inputs.reflectance_calculations_fileName = ['hyperspectral_reflectance_calculations_', char(datetime("today")),'.mat'];





% ------------------------------------------------------
% ----- Define Radiative Transfer Model Parameters -----
% ------------------------------------------------------

% Define the parameters of the INP file


GN_inputs = create_uvSpec_DISORT_inputs_for_EMIT(GN_inputs,emit);










% ----- ISSUE A WARNING! SETTINGS SHOULD BE CHECKED -----

warning([newline, 'Check GN_inputs structure to make sure the settings reflect the situation you wish to model!', newline]);

end
