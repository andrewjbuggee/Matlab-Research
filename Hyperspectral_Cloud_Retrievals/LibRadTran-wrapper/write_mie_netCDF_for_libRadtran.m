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

elseif strcmp(computer_name,'andrewbuggee')==true

    mie_folder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Mie_netCDF_files/'];

elseif strcmp(computer_name,'curc')==true


    mie_folder = '/scratch/alpine/anbu8374/Mie_Calculations/';

end


if ~exist(mie_folder, 'dir')

    mkdir(mie_folder)
end






%%
% ------------------------------------------------------------
% ---------------------- WRITE MIE FILE -----------------------
% ------------------------------------------------------------

wavelengths = 300:5:305;       % HySICS range
radius_groups = [1,10;...
    11,20;...
    21,30;...
    31,40;...
    41,50];

num_wavelengths = length(wavelengths);
num_radius_groups = size(radius_groups, 1);


% Create comments for each line
comments = {'# Mie code to use', '# refractive index to use', '# specify effective radius grid (microns)',...
    '# Specify size distribution and distribution width','# Define wavelength boundaries (nanometers)',...
    '# number of Stokes parameters - just keep total intesnity term', '# maximum number of scattering angles',...
    '# Number of Legrendre terms to be computed (moments of the phase function)',...
    '# multiplicative factor with r_eff which sets distribution upper limit', '# define output variables','# error file length'};

input_filename = cell(num_wavelengths*num_radius_groups, 1);
output_filename = cell(num_wavelengths*num_radius_groups, 1);

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

        % Create the water cloud file
        fileID = fopen([mie_folder, input_filename{idx}], 'w');





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

        % if wavelength has only a single entry, then this is a monochromatic
        fprintf(fileID, '%s          %s \n\n','n_r_max 8 ', comments{9});





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

    runMIE(mie_folder,input_filename{nn},output_filename{nn}, computer_name);

end
toc



%% Stitch the netCDF files into one

