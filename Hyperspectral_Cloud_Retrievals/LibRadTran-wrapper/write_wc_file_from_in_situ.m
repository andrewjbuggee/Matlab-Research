%% This function will write a .DAT water cloud file for LibRadTran

% This function will read and interpolate precomputed Mie calculations for
% water droplets of varrying radii.

% INPUTS:
%   (1) re - effective droplet radius (microns) - this is either a single
%   value, a vector, or a matrix. A single value for re tells the function
%   to create a cloud with a single layer containing a constant droplet
%   radius value. A vector tells the function to create a single wc file
%   with a droplet profile. The length of the vector is equal to the number
%   of layers modeled. A matrix tells the function to create multiple wc
%   files, where the number of columns is equal to the number of wc files
%   created. The number of rows is equal to the number of layers modeled.
%   To create many water cloud files at once that model a homogenous cloud,
%   simply set the column vectors of re to be identical values.
%   ***IMPORTANT*** the re values must start at the cloud bottom with the
%   first value (or row). The last value (or row) is the droplet size at
%   cloud top.



%   (3) z - altitude above sea level (kilometers) - this is a
%   vector with two values: [z_cloudTop, z_cloudBottom]. LibRadTran
%   constructs a cloud by creating layers, where each layer is homogenous
%   untill the next layer is defined. z_cloudTop defines where the cloud
%   ends; this is where the LWC should go to zero. This function will
%   compute a z vector equal in length to that of re using z_topBottom and
%   the geometric thickness H. If re is a matrix, the function expects
%   z_topBottom to be a matrix, where each column is a new wc file. If re
%   is a matrix, and z_topBottom is a vector, then this will be used for
%   every wc file.


%   (5) lambda - wavelength that defines the cloud optical depth
%   (nanometers) - If creating a single wc file, lambda is a single value. If
%   creating multiple wc files, lambda is a vector equal in length to the
%   number of columns of re. If re is a matrix and lambda is a single
%   value, this value will be used for each wc file created.

%   (6) distribution_str - a string telling the code which droplet size
%   distribution to use  - One can chose from two options:
%       (a) 'mono' - monodispersed distribution
%       (b) 'gamma' - gamma droplet distribution.
%       *** IMPORTANT *** For now, this function will NOT
%       use precomputed mie calculations using a gamma droplet
%       distribution. The values returned by LibRadTran appear erroneously
%       high. Instead, if one wishes to use a gamma droplet distribution,
%       the homogenous pre-computed mie table will be used to retrieve mie
%       properties, and then this code will integrate those values over the
%       size distribution.

%   (6) distribution_var - the variance of the size distribution, if
%   applicable. If one is modelling a homogenous cloud, this input will be
%   ignored, and the value can be anything.

%   (7) vert_homogeneity_str - a string telling the code if the cloud is to be
%   modeled as vertically homogeneous. If vertically homogenous, the code will
%   assume every single r_e value represents a single cloud with a constant
%   droplet radius. If vertically non-homogenous, each column of re is assumed
%   to be a single cloud with a droplet profile.
%       (a) 'vert-homogeneous' - homogenous cloud where the entire cloud can be
%       modeled as a single layer with constant properties. If this option
%       is chosen, the code will expect a vector for re, where each entry
%       represents a different cloud.
%       (b) 'vert-non-homogeneous' - a non-homogeneous cloud implies a cloud
%       with multiple layers where the properties vary. If this option is
%       chosen, the code expects a single column vector for re, or a
%       matrix, where each column vector represents a cloud


%   (8) parameterization_str - a string telling the code which
%   parameterization to use when using optical depth and droplet radius to
%   compute the liquid water content. There are 2 options
%       (a) 'mie' - this option uses a pre-computed mie table and
%       interpolates to find values for the extinction efficiency at each
%       wavelength and each droplet radius
%       (b) '2Limit' - this option assumes we are close to the extinction
%       paradox limit and sets the extinction efficiency at a constant
%       value of 2

%   (9) ind_var - the independent variable used to define the effective
%  radius profile. Options are
%           (1) 'optical_depth'
%           (2) 'altitude'

%   (10) compute_weighting_functions - a true or false flag that tells the
%   compute to compute weighting functions by creating N wc files where N
%   is equal to the number of cloud layers. The first file has all N
%   layers, and each successive file removes one cloud layer from the
%   bottom

%   (11) computer_name - the computer this code is running on

%   (12) index - this is the unique identifier that ensures files are not
%   written over one another. If one file is created, the index will be 1.
%   If many fiels are written in a for loop, each file will be tagged with
%   the number in the loop.

%   (13) num_re_parameters - There can be 1 free droplet size parameter for
%   homogeneous clouds, or 2/3 for inhomogeneous clouds. For inhomogeneous,
%   there are either 2 free parameters (re_top and re_bot)
%   or 3 (re_top, re_middle, and re_bot)


%   (14) wc_folder_path - Used to define where the water cloud wc files are
%   store

%   (15) mie_folder_path - Used to define where the mie calculations are
%   stored


% OUTPUTS:
%   (1) .Dat file saved in the libRadTran folder:
%   /.../libRadtran-2.0.4/data/wc

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

%%

function [fileName] = write_wc_file_from_in_situ(re, lwc, z,...
    inSitu_campaign, inSitu_date, inSitu_time, compute_weighting_functions,...
    computer_name, index, wc_folder_path)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 10 inputs

if nargin~=10
    error([newline,'Should be 10 inputs: droplet effective radius, LWC, altitude,',...
        ' in-situ campaign name, date of data recording, time of data recording',...
        ', a flag telling the code to compute weighting functions, the computer name,',...
        'a unique file index,'...
        ' and the location where wc files will be saved.', newline])
end

% Check to make sure re is the same length as the altitude vector and the
% LWC vector

if length(re) ~= length(lwc) || length(re) ~= length(z)

    error([newline,'The Inputs re, lwc, and z must be the same length', newline])

end




% ----- Check to see if there are any NaNs in the radius vector -----

if any(isnan(re))==true

    error([newline, 'The effective radius has atleast one NaN value.', newline])
end







%% Write the WC file


% How many layers to model in the cloud?
nLayers = length(re) + 1;             % Number of altitude levels we need to define a cloud


% Define the filename

fileName = ['WC_in-situ_', inSitu_campaign,'_', char(inSitu_date),'_', num2str(inSitu_time),...
    '-UTC_nn',num2str(index), '.DAT'];


% Determine if the plane is ascending or descending
if (z(1)-z(2))>0
    % the plane is ascending
    % the first few values are at cloud top and there is no change required

elseif (z(1)-z(2))<0
    % the plane is descending
    % the first few values are at cloud bottom and we need to flip re, lwc,
    % and z
    re = flipud(re);
    lwc = flipud(lwc);
    z = flipud(z);

end

% estimate mean dz
dz = -mean(diff(z));   % km


% ------------------------------------------------------------
% ----------------- WE NEED TO APPEND ZEROS ------------------
% ------------------------------------------------------------

% Wherever the cloud is, there needs to be zeros at the cloud top altitude,
% and below the cloud bottom altitude. This information tells LibRadTran
% where the boundaries of the cloud are

% both the effective radius and the LWC need zeros on either boundary,
% unless if the cloud is at the surface


if z(end)==0
    % If true, then the cloud starts at the surface and we only append
    % zeros above the cloud
    re_2write = [re; 0];
    lwc_2write = [lwc; 0];
    z_2write = [z(1) + dz; z];

else
    % In this case, we need zeros below the cloud bottom, and at cloud
    % top
    z_2write = [z(1) + dz; z; 0];                 % create a value at the surface where the cloud parameters go to zero
    re_2write = [0; re; 0];
    lwc_2write = [0; lwc; 0];

end





% ------------------------------------------------------------
% ---------------------- WRITE WC FILE -----------------------
% ------------------------------------------------------------


if compute_weighting_functions==true

    % to compute the weighting functions using the methods of Platnick
    % (2000), we need to compute the reflectance at cloud top (or TOA)
    % by incrementally adding a layer at cloud bottom. Using the z,re,
    % and lwc above, we will iterate through each layer and create n
    % files where n is equal to the number of cloud layers. Each
    % iteration will remove a layer from cloud bottom.

    % compute the optical depth at each layer
    %tau_layer = lwc .* ext_bulk_coeff_per_LWC .* dz_km;       % optical thickness of each homogeneous layer
    % we need a filename for each layer
    filename = cell(nLayers-1, 1);

    for LL = 1:nLayers-1

        % define the filename
        filename{LL} = [fileName{1}(1:end-7), 'layers1-', num2str(nLayers-LL),...
            '.DAT'];
        % Create the water cloud file
        fileID = fopen([wc_folder_path,filename{LL}], 'w');

        % fprintf writes lines in our text file from top to botom
        % wc.DAT files are written with the higher altitudes at the top, and the
        % surface at the bottom

        % to write column vectors in a text file, we have to store them as row
        % vectors

        toWrite = [flipud(z_2write([1, LL+1:end]))'; flipud(lwc_2write([1, LL+1:end]))';...
            flipud(re_2write([1, LL+1:end]))'];

        % Create the opening comment lines of the WC.DAT file

        fprintf(fileID, '%s %10s %7s %8s \n','#','z','LWC','R_eff');
        fprintf(fileID, '%s %10s %7s %8s \n','#','(km)','(g/m^3)','(micron)');

        % Write in the data
        fprintf(fileID,'%12.3f %7.4f %8.3f \n', toWrite);
        fclose(fileID);


    end

    % redefine the output file names
    fileName = filename;


else



    % Create the water cloud file
    fileID = fopen([wc_folder_path,fileName], 'w');

    % fprintf writes lines in our text file from top to botom
    % wc.DAT files are written with the higher altitudes at the top, and the
    % surface at the bottom

    % to write column vectors in a text file, we have to store them as row
    % vectors

    toWrite = [z_2write'; lwc_2write'; re_2write'];

    % Create the opening comment lines of the WC.DAT file

    fprintf(fileID, '%s %10s %7s %8s \n','#','z','LWC','R_eff');
    fprintf(fileID, '%s %10s %7s %8s \n','#','(km)','(g/m^3)','(micron)');

    % Write in the data
    fprintf(fileID,'%12.3f %7.4f %8.3f \n', toWrite);
    fclose(fileID);


end



end





