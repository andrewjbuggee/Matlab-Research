%% Wrapper for Creating Water Cloud files


% By Andrew John Buggee
%%

% Define the folder for where this INP file should be saved
folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/testing_UVSPEC/'];


% Define override of total precipitable water. This will force the total
% column of water vapor to be whatever value you define.
% If you don't wish to change the default, define the variable with nan
total_h2O_column = 10;        % mm of precipitable water


% define the wavelength
wavelength = 555;              % nm - monochromatic wavelength calcualtion

% define the solar zenith angle
sza = 35;                       % degree

% define the viewing zenith angle
vza = 15;                        % degree

% define the atmospheric data file
atm_file = 'afglus.dat';

% Set cell array to hold file names
inputName = cell(1,length(wavelength));
outputName = cell(1,length(wavelength));

% Define empty array for the optical depth
tau_c = zeros(1,length(wavelength));

% define a droplet profile
% define the geometric location of the cloud top and cloud bottom
z_topBottom = [1.5,1];          % km above surface 
z = linspace(z_topBottom(2), z_topBottom(1),25);        % km

re = 10*z - 5;              % microns - droplet profile linear with geometric depth



for nn = 1:length(wavelength)



    % -----------------------------------
    % ---- Write a Water Cloud file! ----
    % -----------------------------------

    % if you wish to have a cloud, define the flag variable below as true
    yesCloud = true;


    % define the total cloud optical depth
    % Round the value to 2 significant digits
    % If you don't, the random number generator will define an optical
    % depth with infinite precision...the digits will continue on for too
    % long, causing issues with dividing the optical depth into layers.
    tau_c(nn) = round(37*rand(1) + 3,2);

%     % define the droplet size
%     % Round the value to 2 significant digits
%     % If you don't, the random number generator will define an effective
%     % radius with infinite precision...
%     re(nn) = round(24*rand(1) + 1,2);                        % microns


    % define the type of droplet distribution
    distribution_str = 'gamma';
    % define the variance of the droplet distribution
    distribution_var = 7;
    % define whether this is a vertically homogenous cloud or not
    vert_homogeneous_str = 'vert-non-homogeneous';
    % define how liquid water content will be computed
    parameterization_str = 'mie';



    wc_filename = write_wc_file(re,tau_c(nn),z_topBottom, wavelength(nn), distribution_str,...
        distribution_var, vert_homogeneous_str, parameterization_str);

    % Define the percentage of cloud cover
    percent_cloud_cover = 1;

    % Define the parameterization scheme used to comptue the optical quantities
    wc_parameterization = 'mie interpolate';
    % --------------------------------------------------------------
    % --------------------------------------------------------------





    % define the name of the input file
    if yesCloud==true
        inputName{nn} = [num2str(wavelength(nn)),'nm_withCloudLayer_',num2str(sza),'sza_',num2str(vza),...
            'vza_',atm_file(1:end-4),'.INP'];
    else
        inputName{nn} = [num2str(wavelength(nn)),'nm_',num2str(sza),'sza_',num2str(vza),...
            'vza_',atm_file(1:end-4),'.INP'];
    end



    % Define the .OUT file name
    outputName{nn} = ['OUTPUT_',inputName{nn}(1:end-4)];




    % ---- ****************** ------
    % ---- Write the INP File ------
    % ---- ****************** ------

    % Open the old file for writing
    fileID = fopen([folderName,inputName{nn}], 'w');

    % Define which RTE solver to use
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'rte_solver','disort',' ', '# Radiative transfer equation solver');


    % Define the number of streams to keep track of when solving the equation
    % of radiative transfer
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'number_of_streams','6',' ', '# Number of streams');


    % Define the location and filename of the atmopsheric profile to use
    % ------------------------------------------------
    formatSpec = '%s %5s %s \n';
    fprintf(fileID, formatSpec,['atmosphere_file ','../data/atmmod/',atm_file],' ', '# Location of atmospheric profile');

    % Define the location and filename of the extraterrestrial solar source
    % ---------------------------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'source solar','../data/solar_flux/kurudz_1.0nm.dat', ' ', '# Bounds between 250 and 10000 nm');


    % Define the total precipitable water
    % ------------------------------------------------
    if isnan(total_h2O_column)==false
        formatSpec = '%s %s %f %s %5s %s \n';
        fprintf(fileID, formatSpec,'mol_modify','H2O', total_h2O_column,' MM', ' ', '# Total Precipitable Water');
    end


    % Define the surface albedo
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n\n';
    fprintf(fileID, formatSpec,'albedo','0', ' ', '# Surface albedo of the ocean');


    if yesCloud==true

        % Define the water cloud file
        % ------------------------------------------------
        formatSpec = '%s %s %5s %s \n';
        fprintf(fileID, formatSpec,'wc_file 1D', ['../data/wc/',wc_filename{1}], ' ', '# Location of water cloud file');


        % Define the percentage of horizontal cloud cover
        % This is a number between 0 and 1
        % ------------------------------------------------
        formatSpec = '%s %f %5s %s \n';
        fprintf(fileID, formatSpec,'cloudcover wc', percent_cloud_cover, ' ', '# Cloud cover percentage');


        % Define the technique or parameterization used to convert liquid cloud
        % properties of r_eff and LWC to optical depth
        % ----------------------------------------------------------------------
        formatSpec = '%s %s %5s %s \n\n';
        fprintf(fileID, formatSpec,'wc_properties', wc_parameterization, ' ', '# optical properties parameterization technique');

    end

    % Define the wavelengths for which the equation of radiative transfer will
    % be solve
    % -------------------------------------------------------------------------
    formatSpec = '%s %f %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'wavelength', wavelength(nn), wavelength(nn), ' ', '# Wavelength range');


    % Define the sensor altitude
    % ------------------------------------------------
    formatSpec = '%s %s %5s %s \n';
    fprintf(fileID, formatSpec,'zout', 'toa', ' ', '# Sensor Altitude');

    % Define the solar zenith angle
    % ------------------------------------------------
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'sza', sza, ' ', '# Solar zenith angle');

    % Define the solar azimuth angle
    % ------------------------------------------------
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'phi0', 0.0, ' ', '# Solar azimuth angle');

    % Define the cosine of the zenith viewing angle
    % ------------------------------------------------
    formatSpec = '%s %f %5s %s \n';
    fprintf(fileID, formatSpec,'umu', cosd(vza), ' ', '# Cosine of the zenith viewing angle');

    % Define the azimuth viewing angle
    % ------------------------------------------------
    formatSpec = '%s %f %5s %s \n\n';
    fprintf(fileID, formatSpec,'phi', 0.0, ' ', '# Azimuthal viewing angle');


    % Set the error message to quiet of verbose
    % ------------------------------------------------
    formatSpec = '%s';
    fprintf(fileID, formatSpec,'verbose');


    % Close the file!
    fclose(fileID);
    % ----------------------------------------------------
    % ----------------------------------------------------


end


%%  Run the created INP files!


for nn = 1:length(wavelength)
    runUVSPEC(folderName,inputName{nn},outputName{nn});
end


