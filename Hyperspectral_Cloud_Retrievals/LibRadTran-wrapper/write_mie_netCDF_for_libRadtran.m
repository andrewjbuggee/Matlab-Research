%% Write mie files that generate a netCDF output. Loop over wavelength and we will stitch
% together each cdf file at the end


% Andrew John Buggee

clear variables

%% Find folders

computer_name = whatComputer;

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(computer_name,'anbu8374')==true

    mie_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_netCDF_files/';

    master_folder = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_netCDF_files/master_netCDF_folder/';

elseif strcmp(computer_name,'andrewbuggee')==true

    mie_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Mie_netCDF_files/'];

    master_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Mie_netCDF_files/master_netCDF_folder/'];

elseif strcmp(computer_name,'curc')==true


    mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/Mie_netCDF_files/';

    master_folder = '/scratch/alpine/anbu8374/Mie_Calculations/master_netCDF_folder/';

end


if ~exist(mie_folder, 'dir')

    mkdir(mie_folder)
end

if ~exist(master_folder, 'dir')

    mkdir(master_folder)
end






%%
% ------------------------------------------------------------
% ---------------------- WRITE MIE FILE -----------------------
% ------------------------------------------------------------

wavelengths = 300:5:305;       % HySICS range

radius_range = [1, 50];
num_radius_groups = 10;
radius_groups = zeros(num_radius_groups, 2);
for nn = 1:num_radius_groups
    if nn == 1
        radius_groups(nn,:) = [1, nn*(radius_range(end)/num_radius_groups)];
    else
        radius_groups(nn,:) = [radius_groups(nn-1,2)+1, nn*(radius_range(end)/num_radius_groups)];
    end
end


% wavelengths = 300:5:305;       % HySICS range
% radius_groups = [1,2];

num_wavelengths = length(wavelengths);


% Create comments for each line
comments = {'# Mie code to use', '# refractive index to use', '# specify effective radius grid (microns)',...
    '# Specify size distribution and distribution width','# Define wavelength boundaries (nanometers)',...
    '# number of Stokes parameters - just keep total intesnity term', '# maximum number of scattering angles',...
    '# Number of Legrendre terms to be computed (moments of the phase function)',...
    '# multiplicative factor with r_eff which sets distribution upper limit', '# define output variables','# error file length'};

input_filename = cell(num_wavelengths*num_radius_groups, 1);
output_filename = cell(num_wavelengths*num_radius_groups, 1);
foldername = cell(num_wavelengths*num_radius_groups, 1);

idx = 0;

changing_variables = zeros(num_radius_groups*num_wavelengths, 3);

for rr = 1:num_radius_groups

    for ww = 1:num_wavelengths

        idx = idx+1;

        % store the variable configurations for each iteration
        changing_variables(idx,:) = [radius_groups(rr,:), wavelengths(ww)];

        % --------------------------------------------
        % *** create the input and output filename ***
        % --------------------------------------------

        input_filename{idx} = ['Mie_calc_water_GammaDist_HySICS_',num2str(wavelengths(ww)),...
            'nm_re_', num2str(radius_groups(rr, 1)), '-', num2str(radius_groups(rr,2)),'microns.INP'];
        output_filename{idx} = ['OUTPUT_', input_filename{idx}(1:end-4)];


        % --------------------------------------------
        % *** create a unique folder for each file ***
        % --------------------------------------------
        % This has to be done because otherwise libRadtran writes over each
        % netCDF file created
        foldername{idx} = [mie_folder,input_filename{idx}(1:end-4), '/'];

        if ~exist(foldername{idx}, 'dir')

            mkdir(foldername{idx})

        end



        % Create the water cloud file
        fileID = fopen([foldername{idx}, input_filename{idx}], 'w');





        % fprintf writes lines in our text file from top to botom
        % .INP files for mie calculations always require the same inputs

        % to write column vectors in a text file, we have to store them as row
        % vectors

        % ----------------------------------
        % Define the mie program code to use
        % ----------------------------------

        fprintf(fileID, '%s          %s \n','mie_program MIEV0 ', comments{1});



        % ----------------------------------
        % Define the Index of Refraction!
        % ----------------------------------

        % check to see if the index of refraction is a string or a number
        fprintf(fileID, '%s          %s \n','refrac water ', comments{2});



        % ---------------------------------------------------------------------
        % Write in the value for the modal radius. Check to see if its a vector
        % ---------------------------------------------------------------------


        fprintf(fileID,'%s %f %f %s          %s \n', 'r_eff', radius_groups(rr,1), radius_groups(rr,2), '1', comments{3});







        % ----------------------------------------------------------------
        % Write in the value for the droplet distribution, if its not mono
        % ----------------------------------------------------------------



        fprintf(fileID,'%s         %s \n', 'distribution gamma 7', comments{4});




        % ---------------------------------------
        % ----- define the wavelength range -----
        % ---------------------------------------

        % if wavelength has only a single entry, then this is a monochromatic
        fprintf(fileID,'%s  %f %f          %s \n\n', 'wavelength', wavelengths(ww), wavelengths(ww), comments{5});



        % ---------------------------------------
        % ----- define stokes parameters -----
        % ---------------------------------------

        % if wavelength has only a single entry, then this is a monochromatic
        fprintf(fileID, '%s          %s \n','nstokes 1 ', comments{6});


        % ---------------------------------------
        % ----- maximum number of scattering angles -----
        % ---------------------------------------

        % if wavelength has only a single entry, then this is a monochromatic
        fprintf(fileID, '%s          %s \n','nthetamax 1000 ', comments{7});


        % ---------------------------------------
        % ----- define the number of legendre polynomials to compute -----
        % ---------------------------------------

        % if wavelength has only a single entry, then this is a monochromatic
        fprintf(fileID, '%s          %s \n','nmom 10000 ', comments{8});


        % ---------------------------------------
        % ----- define the multiplicative factor for the upper range of the sixe distribution -----
        % ---------------------------------------

        % How many radii do you wish to calculated when defining the
        % droplet distribution? 8 takes a long time!
        fprintf(fileID, '%s          %s \n\n','n_r_max 5 ', comments{9});





        % ---------------------------
        % Define the output variables
        % ---------------------------

        % But first make a comment
        fprintf(fileID,'\n%s \n', comments{10});
        fprintf(fileID,'%s \n', 'output_user netcdf');





        % -----------------------
        % Print the error message
        % -----------------------
        fprintf(fileID, '%s          %s', 'verbose', comments{11});


        % Close the file!
        fclose(fileID);



    end


end



%% Run the mie files!

tic
parfor nn = 1:length(input_filename)

    disp(['nn = ', num2str(nn), newline])

    runMIE(foldername{nn},input_filename{nn},output_filename{nn}, computer_name);

end
toc



%% step through each folder, rename the netCDF file, and move it into a
% master folder


for nn = 1:length(foldername)

    % Find all files that end with .cdf
    contents = dir([foldername{nn}, '*.cdf']);

    if isempty(contents)==false


        % Rename and move the netCDF file
        old_folder_file = [foldername{nn}, contents.name];
        new_folder_file = [master_folder, input_filename{nn}(1:end-4),'.mie.cdf'];

        cmd = ['mv ', old_folder_file, ' ', new_folder_file];
        system(cmd);


    end


end





%% Take 3 - I think I'll have to do it myslef

contents_master = dir([master_folder, '*.cdf']);

% info1 = ncinfo([contents_master(1).folder,'/', contents_master(1).name]);
% info2 = ncinfo([contents_master(2).folder,'/', contents_master(2).name]);

% Concatenate the vairables from each file
% Each file computes the mie variables for a unique r_eff and wavelength

% this variable is the spread of wavelengths according to 
wavelen = wavelengths;
reff = (1:50)';

% theta is a 4-D matrix [angles, stokes parameters, effective radii,wavelengths]
% I've hard coded my mie files to compute 1000 scattering angles, the first
% stokes parameter, 50 effective radii over some set of wavelengths
theta = zeros(1000, 1, 50, length(wavelengths));

% Number of scattering angles is a 3-D matrix
% [stokes parameters, effective radii,wavelengths]
ntheta = zeros(1, 50, length(wavelengths));

% Scattering Phase matrix
% a 4-D matrix [angles, stokes parameters, effective radii, wavelengths]
% I've hard coded my mie files to compute 1000 scattering angles, the first
% stokes parameter, 50 effective radii over some set of wavelengths
phase = zeros(1000, 1, 50, length(wavelengths));

% Legendre Polynomials
% a 4-D matrix [moments, stokes parameters, effective radii, wavelengths]
% I've hard coded my mie files to compute 10000 polynomials, the first
% stokes parameter, 50 effective radii over some set of wavelengths
pmom = zeros(10000, 1, 50, length(wavelengths));

% The number of Legendre polynomials
% a 3-D matrix [stokes parameters, effective radii, wavelengths]
% Only computing the first stokes parameter so get rid of the first
% dimension
% The maximum possible number of Legendre polynomials my files will store
% is 10,000. But they may not all compute this many!
nmom = zeros(1, 50, length(wavelengths));

% Extinction efficiency
% a 2-D matrix [effective radii, wavelengths]
ext = zeros(50, length(wavelengths));

% Single Scattering Albedo
% a 2-D matrix [effective radii, wavelengths]
ssa = zeros(50, length(wavelengths));

% Asymmetry Parameter
% a 2-D matrix [effective radii, wavelengths]
gg = zeros(50, length(wavelengths));

% Real part of the refractive index
% varies with number of wavelengths
% row vector 
refre = zeros(length(wavelengths), 1);

% Real part of the refractive index
% varies with number of wavelengths
% row vector 
refim = zeros(length(wavelengths), 1);

% Density of liquid water (g/cm^3)
rho = 1;


for nn = 1:length(contents_master)

    info = ncinfo([contents_master(nn).folder,'/', contents_master(nn).name]);

    % This is just one value per file - raw value is in microns
    wavelen_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'wavelen') * 1e3;  % nanometers

    % find the wavelength indexes
    wl_idx = find(wavelen_temp==wavelen);

    % read these values to know where the variables below need to be
    % indexed
    reff_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'reff');

    % find the effective radii indexes
    % Because the effective radii range from 1-50 microns, the indexes are
    % simply the effective radii computed for the file in question
    re_idx = reff_temp;

    % % find the effective radii indexes
    % clear re_idx
    % if length(reff_temp)>1
    % 
    %     for rr = 1:length(reff_temp)
    % 
    %         re_idx(rr) = find(reff_temp(rr)==reff);
    %     end
    % 
    % elseif isscalar(reff_temp)
    %     re_idx = find(reff_temp==reff);
    % end


    

    % Theta never changes ?
    theta_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'theta');

    % ntheta deffinitely donest change, this is hard coded
    ntheta_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'ntheta');


    phase_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'phase');
    pmom_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'pmom');
    nmom_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'nmom');
    ext_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'ext');
    ssa_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'ssa');
    refre_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'refre');
    refim_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'refim');
    rho_temp = ncread([contents_master(nn).folder,'/', contents_master(nn).name], 'rho');



    % --- Concatenate! ---

    % Theta is a four dimensional variable of scattering angles
    

    % concatenate ntheta
    ntheta(1, re_idx, wl_idx) = ntheta_temp;

    % concatenate nmom
    nmom(1, re_idx, wl_idx) = nmom_temp;






    
end




% % Write the schema for a new .CDF file using the first file schema
% ncwriteschema('test_merge.cdf',info1)
% 
% % add the schema from the second file
% ncwriteschema('test_merge.cdf',info2)

%%

% The code creates 3 NetCDF Files and combine them.
numFiles = 3;
dimSize = 10;
% Step 1: Create Multiple NetCDF Files
for i = 1:numFiles
    filename = sprintf('test_file_%d.nc', i);
    nccreate(filename, 'data', 'Dimensions', {'x', dimSize, 'y', dimSize})
    data = rand(dimSize, dimSize);
    ncwrite(filename, 'data', data);
end
% Step 2: Combine NetCDF Files into One
outputFile = 'combined_file.nc';
% 'data' is the variable name in the NetCDF file
nccreate(outputFile, 'data', 'Dimensions', {'x', dimSize, 'y', dimSize, 'file', numFiles});
for i = 1:numFiles
    filename = sprintf('test_file_%d.nc', i);
    data = ncread(filename, 'data');
    ncwrite(outputFile, 'data', data, [1, 1, i]);
end

