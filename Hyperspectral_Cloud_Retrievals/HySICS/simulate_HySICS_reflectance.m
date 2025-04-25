%% Create simulated HySICS measurements
% Include Gaussian random noise


% By Andrew John Buggee


clear variables


%% Define the path location where INP files will be stored, and where Reflectances will be stored

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Simulated_spectra/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/';




elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Simulated_spectra/'];


    % Define the folder path where all .INP files will be saved
    folderpath_inp = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/HySICS/'];

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/data/'];






elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the folder path where .mat files of relfectance will be stored
    folderpath_reflectance = '/scratch/alpine/anbu8374/Thermodynamic_phase/';



    % Define the folder path where all .INP files will be saved
    folderpath_inp = '/scratch/alpine/anbu8374/Thermodynamic_phase/';

    % Define the libRadtran data files path. All paths must be absolute in
    % the INP files for libRadtran
    libRadtran_data_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/';


end


% If the folder path doesn't exit, create a new directory
if ~exist(folderpath_inp, 'dir')

    mkdir(folderpath_inp)

end


% If the folder path doesn't exit, create a new directory
if ~exist(folderpath_reflectance, 'dir')

    mkdir(folderpath_reflectance)

end



%% ---- First, let's simulate water clouds ----


% Define the parameters of the INP file
clear inputs



% Define the RTE Solver
inputs.RT.rte_solver = 'disort';


% Define the number of streams to use in your radiative transfer model
inputs.RT.num_streams = 16;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------- Define the source file and resolution -------------------

% source_file = 'hybrid_reference_spectrum_p1nm_resolution_c2022-11-30_with_unc.dat';
% source_file_resolution = 0.025;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
% inputs.RT.source_file = 'hybrid_reference_spectrum_1nm_resolution_c2022-11-30_with_unc.dat';
% inputs.RT.source_file_resolution = 0.1;         % nm

% these data have 0.1nm sampling resolution, despite what the file name
% suggests
inputs.RT.source_file = 'kurudz_1.0nm.dat';
inputs.RT.source_file_resolution = 1;         % nm

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% ---------------------- DEFINE THE WAVELENGTHS! -------------------------
% ------------------------------------------------------------------------
% define the wavelength range. If monochromatic, enter the same number
% twice


% ----------------- Simulating HySICS spectral channels ------------------
% number of channels = 460
inputs.bands2run = (1:10:460)';



% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Create the spectral response functions ------------------
% ------------------------------------------------------------------------

% ------------------------------------
% modeling the HySICS instrument...
% ------------------------------------
% Define the HySICS spectral response functions
spec_response = create_HySICS_specResponse(inputs.bands2run, inputs.RT.source_file_resolution);

% now define the wavelength range of each spectral channel
inputs.RT.wavelengths2run = zeros(length(inputs.bands2run), 2);

for ww = 1:length(inputs.bands2run)
    % The wavelength vector for libRadTran is simply the lower and upper
    % bounds
    inputs.RT.wavelengths2run(ww,:) = [spec_response{ww}(1, 1),...
        spec_response{ww}(end, 1)];

end


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% --- Do you want to use the Nakajima and Tanka radiance correction? -----
inputs.RT.use_nakajima_phaseCorrection = true;
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ----------------- What band model do you want to use? ------------------

% reptran coarse is the default
% if using reptran, provide one of the following: coarse (default), medium
% or fine
inputs.RT.band_parameterization = 'reptran coarse';
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% define the atmospheric data file
inputs.RT.atm_file = 'afglus.dat';

% define the surface albedo
inputs.RT.albedo = 0.04;

% day of the year
inputs.RT.day_of_year = 100;
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% -------------- Do you want a cloud in your model? ----------------------
inputs.RT.yesCloud = true;

% define the cloud geometric depth
inputs.RT.cloud_depth = 500;                % meters

% define the geometric location of the cloud top and cloud bottom
inputs.RT.z_topBottom = [1.5, 1];          % km above surface


% Water Cloud depth
inputs.RT.H = inputs.RT.z_topBottom(1) - inputs.RT.z_topBottom(2);                                % km - geometric thickness of cloud
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% --------------------- Various Cloud modeling inputs --------------------
% ------------------------------------------------------------------------
% Do you want use your custom mie calculation file?
inputs.RT.use_custom_mie_calcs = false;


% define how liquid water content will be computed in the write_wc_file
% function.
inputs.RT.parameterization_str = 'mie';

% define the wavelength used for the optical depth as the 650 nm
% band1 = modisBands(1);
% lambda_forTau = band1(1);            % nm
inputs.RT.lambda_forTau = 500;            % nm
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% -------------------- Cloud optical properties --------------------------
% ------------------------------------------------------------------------

% define whether this is a vertically homogenous cloud or not
inputs.RT.vert_homogeneous_str = 'vert-non-homogeneous';


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % --------- HOMOGENOUS CLOUD -------------



    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.dist_var = 7;              % distribution variance

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    inputs.RT.re = 10.79;      % microns
    inputs.RT.tau_c = 7;


elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------- NON-HOMOGENOUS CLOUD -------------

    % For a cloud with a nonconstant droplet profile
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

    inputs.RT.profile_type = 'adiabatic'; % type of water droplet profile

    inputs.RT.n_layers = 10;                          % number of layers to model within cloud

    inputs.RT.z = linspace(inputs.RT.z_topBottom(1), inputs.RT.z_topBottom(2), inputs.RT.n_layers);        % km - altitude above ground vector

    inputs.RT.indVar = 'altitude';                    % string that tells the code which independent variable we used

    inputs.RT.dist_var = linspace(10,10, inputs.RT.n_layers);              % distribution variance

    inputs.RT.r_top = 12.565;     % microns
    inputs.RT.r_bot = 4.135;        % microns
    inputs.RT.tau_c = 6.424;

    % define the type of droplet distribution
    inputs.RT.distribution_str = 'gamma';

    % define the spread of the droplet distribution
    % *** To match the optical
    %   properties mie table precomputed by libRadtran, use a gamma
    %   distribution alpha parameter of 7 ***
    inputs.RT.dist_var = 7;              % distribution variance


end

% ------------------------------------------------------------------------



% Define the parameterization scheme used to comptue the optical quantities
% within the INP file, i.e. what function call that tells libRadtran how to
% compute scattering and optical quantities
if inputs.RT.use_custom_mie_calcs==false

    inputs.RT.wc_parameterization = 'mie interpolate';
    %inputs.RT.wc_parameterization = 'hu';

else
    %wc_parameterization = '../data/wc/mie/wc.mie_test.cdf interpolate';
    inputs.RT.wc_parameterization = '../data/wc/mie/wc.mie_test2_more_nmom.cdf interpolate';
end

% --------------------------------------------------------------
% --------------------------------------------------------------



% --------------------------------------------------------------
% ----------- Define the Solar and Viewing Gemometry -----------
% --------------------------------------------------------------

% Define the altitude of the sensor
inputs.RT.sensor_altitude = 'toa';          % top-of-atmosphere

% define the solar zenith angle
inputs.RT.sza = 19.5688;           % degree

% Define the solar azimuth measurement between values 0 and 360
% The EMIT solar azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran solar azimuth is defined as 0-360 degrees
% clockwise from due south. So they are separated by 180 degrees. To map
% the EMIT azimuth the the libRadTran azimuth, we need to add 180 modulo
% 360
inputs.RT.phi0 = 113.8140;         % degree

% define the viewing zenith angle
inputs.RT.vza = 8.3134; % values are in degrees;                        % degree

% define the viewing azimuth angle
% The EMIT sensor azimuth angle is defined as 0-360 degrees clockwise from
% due north. The libRadTran sensor azimuth is defined as 0-360 degrees
% clockwise from due North as well. So they are separated by 180 degrees. A
% sensor azimuth angle of 0 means the sensor is in the North, looking
% south. No transformation is needed

inputs.RT.vaz = 70.0849;     % degree
% --------------------------------------------------------------



% --------------------------------------------------------------
% --- Do you want to use the Cox-Munk Ocean Surface Model? -----
inputs.RT.use_coxMunk = true;
inputs.RT.wind_speed = 3;             % m/s
% --------------------------------------------------------------


% ------------------------------------------------------------------------
% --------- Do you want boundary layer aerosols in your model? -----------
inputs.RT.yesAerosols = true;

inputs.RT.aerosol_type = 4;               % 4 = maritime aerosols
inputs.RT.aerosol_opticalDepth = 0.1;     % MODIS algorithm always set to 0.1
% ------------------------------------------------------------------------




% --------------------------------------------------------
% --------- What is column water vapor amount? -----------

% Use a custom H2O profile
inputs.RT.H2O_profile = 'afglus_H2O_none_inside_cloud.dat';


% Using measurements from the AMSR2 instrument, a passive microwave
% radiometer for 17 Jan 2024
inputs.RT.modify_waterVapor = false;

inputs.RT.waterVapor_column = 20;              % mm - milimeters of water condensed in a column
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% ------- Do you want to modify concentration of Carbon dioxide? ---------

% 400 ppm = 1.0019 * 10^23 molecules/cm^2
inputs.RT.modify_CO2 = true;

inputs.RT.CO2_mixing_ratio = 410;       % ppm
% ------------------------------------------------------------------------


% --------------------------------------------------------------
% --- Do you want to uvSpec to compute reflectivity for you? ---
inputs.RT.compute_reflectivity_uvSpec = false;
% --------------------------------------------------------------



% --------------------------------------------------------------
% Do you want to print an error message?
inputs.RT.errMsg = 'verbose';
% --------------------------------------------------------------


%% Write each INP file



% num wavelengths
if size(inputs.RT.wavelengths2run,1)>1 && size(inputs.RT.wavelengths2run,2)>1

    inputs.RT.monochromatic_calc = false;
    num_wl = size(inputs.RT.wavelengths2run,1);

else

    inputs.RT.monochromatic_calc = true;
    num_wl = length(inputs.RT.wavelengths2run);

end


idx = 0;


tic

if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------

    inputFileName = cell(length(inputs.RT.re), length(inputs.RT.tau_c), num_wl);
    outputFileName = cell(length(inputs.RT.re), length(inputs.RT.tau_c), num_wl);


    for rr = 1:length(inputs.RT.re)



        for tc = 1:length(inputs.RT.tau_c)


            idx = idx + 1;
            % -----------------------------------
            % ---- Write a Water Cloud file! ----
            % -----------------------------------

            % ------------------------------------------------------
            % --------------------VERY IMPORTANT ------------------
            % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
            % ------------------------------------------------------
            wc_filename = write_wc_file(inputs.RT.re(rr), inputs.RT.tau_c(tc), inputs.RT.z_topBottom,...
                inputs.RT.lambda_forTau, inputs.RT.distribution_str, inputs.RT.dist_var,...
                inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, idx);

            inputs.RT.wc_filename = wc_filename{1};


            disp(['Iteration: [re, tc] = [', num2str(rr), '/', num2str(length(inputs.RT.re)),', ',...
                num2str(tc), '/', num2str(length(inputs.RT.tau_c)), ']', newline])


            for ww = 1:num_wl



                % ---- Define the wavelengths ----
                inputs.RT.wavelengths = inputs.RT.wavelengths2run(ww,:);


                disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])


                % ------------------------------------------------
                % ---- Define the input and output filenames! ----
                % ------------------------------------------------
                % input_names need a unique identifier. Let's give them the nn value so
                % they can be traced, and are writen over in memory

                inputFileName{rr, tc, ww} = [num2str(inputs.RT.wavelengths(ww, 1)), '-', num2str(inputs.RT.wavelengths(ww, 2)), 'nm_',...
                    inputs.RT.atm_file(1:end-4),'.INP'];



                outputFileName{rr, tc, ww} = ['OUTPUT_',inputFileName{rr,tc,ww}(1:end-4)];


                % ------------------ Write the INP File --------------------
                write_INP_file(folderpath_inp, libRadtran_data_path, inputFileName{rr, tc, ww}, inputs);




            end



        end

    end





elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------

    inputFileName = cell(length(inputs.RT.r_top), length(inputs.RT.r_bot), length(inputs.RT.tau_c), num_wl);
    outputFileName = cell(length(inputs.RT.r_top), length(inputs.RT.r_bot), length(inputs.RT.tau_c), num_wl);



    % create a legend string
    lgnd_str = cell(1, length(inputs.RT.waterVapor_column));

    % store the reflectances
    Refl_model = zeros(num_wl, length(inputs.RT.r_top), length(inputs.RT.r_bot),...
        length(inputs.RT.tau_c));


    for rt = 1:length(inputs.RT.r_top)


        for rb = 1:length(inputs.RT.r_bot)


            for tc = 1:length(inputs.RT.tau_c)


                disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])


                idx = idx + 1;
                % -----------------------------------
                % ---- Write a Water Cloud file! ----
                % -----------------------------------

                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
                % ------------------------------------------------------

                % -----------------------------------
                % ---- Write a Water Cloud file! ----
                % -----------------------------------
                % most uncertainties for the modis optical retrieval are between 2
                % and 10 percent. So lets round off all re values to the 1000th decimal
                % place

                re = create_droplet_profile2([inputs.RT.r_top(rt), inputs.RT.r_bot(rb)],...
                    inputs.RT.z, inputs.RT.indVar, inputs.RT.profile_type);     % microns - effective radius vector


                % ------------------------------------------------------
                % --------------------VERY IMPORTANT ------------------
                % ADD THE LOOP VARIABLE TO THE WC NAME TO MAKE IT UNIQUE
                % ------------------------------------------------------
                wc_filename = write_wc_file(re, inputs.RT.tau_c(tc), inputs.RT.z_topBottom,...
                    inputs.RT.lambda_forTau, inputs.RT.distribution_str, inputs.RT.dist_var,...
                    inputs.RT.vert_homogeneous_str, inputs.RT.parameterization_str, rt*rb);

                inputs.RT.wc_filename = wc_filename{1};



                for ww = 1:num_wl


                    % ---- Define the wavelengths ----
                    inputs.RT.wavelengths = inputs.RT.wavelengths2run(ww,:);



                    disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])



                    % ------------------------------------------------
                    % ---- Define the input and output filenames! ----
                    % ------------------------------------------------
                    % input_names need a unique identifier. Let's give them the nn value so
                    % they can be traced, and are writen over in memory

                    inputFileName{rt, rb, tc, ww} = [num2str(inputs.RT.wavelengths2run(ww)), 'nm_',...
                        inputs.RT.atm_file(1:end-4),'.INP'];



                    outputFileName{rt, rb, tc, ww} = ['OUTPUT_',inputFileName{rt,rb,tc,ww}(1:end-4)];




                    % ------------------ Write the INP File --------------------
                    write_INP_file(folderpath_inp, libRadtran_data_path, inputFileName{rt, rb, tc, ww}, inputs);





                end



            end

        end

    end


end

toc






%% Calculate Reflectance



tic


if strcmp(inputs.RT.vert_homogeneous_str, 'vert-homogeneous') == true

    % ----------------------------------------
    % --------- HOMOGENOUS CLOUD -------------
    % ----------------------------------------


    % create a legend string
    lgnd_str = cell(1, length(inputs.RT.waterVapor_column));

    % store the reflectances
    Refl_model = zeros(num_wl, length(inputs.RT.re), length(inputs.RT.tau_c));


    for rr = 1:length(inputs.RT.re)



        for tc = 1:length(inputs.RT.tau_c)



            disp(['Iteration: [re, tc] = [', num2str(rr), '/', num2str(length(inputs.RT.re)),', ',...
                num2str(tc), '/', num2str(length(inputs.RT.tau_c)), ']', newline])


            parfor ww = 1:size(inputs.RT.wavelengths2run, 1)
             %for ww = 1:size(inputs.RT.wavelengths2run, 1)


                disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])


                % ----------------------------------------------------
                % --------------- RUN RADIATIVE TRANSFER -------------
                % ----------------------------------------------------


                % compute INP file
                [inputSettings] = runUVSPEC(folderpath_inp, inputFileName{rr, tc, ww}, outputFileName{rr, tc, ww});

                % read .OUT file
                % radiance is in units of mW/nm/m^2/sr
                [ds,~,~] = readUVSPEC(folderpath_inp, outputFileName{rr, tc, ww},inputSettings(2,:),...
                    inputs.RT.compute_reflectivity_uvSpec);

                % Store the Radiance
                %            Rad_model(rr, tc, ww, :) = ds.radiance.value;       % radiance is in units of mW/nm/m^2/sr

                % compute the reflectance
                % Refl_model(ww, rr, tc) = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
                %     spec_response.value(ww, :)');

                [Refl_model(ww, rr, tc), ~] = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));








            end



        end

    end





elseif strcmp(inputs.RT.vert_homogeneous_str, 'vert-non-homogeneous') == true

    % --------------------------------------------
    % --------- NON-HOMOGENOUS CLOUD -------------
    % --------------------------------------------

    % create a legend string
    lgnd_str = cell(1, length(inputs.RT.waterVapor_column));

    % store the reflectances
    Refl_model = zeros(num_wl, length(inputs.RT.r_top), length(inputs.RT.r_bot),...
        length(inputs.RT.tau_c));


    for rt = 1:length(inputs.RT.r_top)


        for rb = 1:length(inputs.RT.r_bot)


            for tc = 1:length(inputs.RT.tau_c)


                disp(['Iteration: [rt, rb, tc] = [', [num2str(rt),', ', num2str(rb), ', ', num2str(tc)], ']...', newline])



                parfor ww = 1:length(inputs.RT.wavelengths2run)
                %for ww = 1:length(inputs.RT.wavelengths2run)



                    disp(['ww = ', num2str(ww),'/',num2str(num_wl), newline])



                    % ----------------------------------------------------
                    % --------------- RUN RADIATIVE TRANSFER -------------
                    % ----------------------------------------------------


                    % compute INP file
                    [inputSettings] = runUVSPEC(folderpath_inp, inputFileName{rt, rb, tc, ww}, outputFileName{rt, rb, tc, ww});

                    % read .OUT file
                    % radiance is in units of mW/nm/m^2/sr
                    [ds,~,~] = readUVSPEC(folderpath_inp, outputFileName{rt, rb, tc, ww}, inputSettings(2,:),...
                        inputs.RT.compute_reflectivity_uvSpec);


                    % compute the reflectance
                    [Refl_model(ww, rt, rb, tc), ~] = reflectanceFunction(inputSettings(2,:), ds, spec_response{ww}(:,2));









                end



            end

        end

    end


end



toc


%% Add Gaussian random noise



%% Plot the results


figure;
if size(inputs.RT.wavelength,1)>1 && size(inputs.RT.wavelength,2)>1

    plot(mean(inputs.RT.wavelength,2), Refl_model, '-', 'linewidth', 5, 'markersize', 27, 'Color', mySavedColors(1, 'fixed'),...
        'linewidth', 3)

else

    plot(inputs.RT.wavelength, Refl_model, '-', 'linewidth', 5, 'markersize', 27, 'Color', mySavedColors(1, 'fixed'),...
        'linewidth', 3)

end

hold on

grid on; grid minor
xlabel('Wavelength (nm)','Interpreter', 'latex')
ylabel('Reflectance (1/sr)','Interpreter', 'latex')
set(gcf, 'Position', [0 0 1000 1000])
legend('Simulated Reflectance',  'Interpreter', 'latex', 'Fontsize', 30', 'location', 'best')
% title(['Simulated Reflectance - liquid water cloud - $r_e = $', num2str(re_2plot), ' $\mu m$, $\tau_c = $',...
%     num2str(tau_2plot)], 'Interpreter', 'latex')
