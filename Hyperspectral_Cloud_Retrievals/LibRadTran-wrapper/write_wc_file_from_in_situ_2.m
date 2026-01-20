%% This function will write a .DAT water cloud file for LibRadTran

% This function uses in-situ measurements of effective raidius and liquid
% water content as a function of altitude above sea level and creates a wc
% file

% INPUTS:
%   (1) re - effective droplet radius (microns) - this is a vector
%   containing the measurements of effective radius as the plane flew
%   through cloud sampling the vertical distribution of droplets 


%   (2) lwc - liquid water content (g/m^3) - this is a
%   vector of the liquid water content sampled along the in-situ flight
%   path


%   (3) z - altitude above sea level (kilometers) - this is a
%   vector of the altitudes at which the effective radius and liquid water
%   content were samples by the in-situ probe

%   (4) alpha - This is the width parameter of the droplet size
%   distribution. We fit gamma distributions according to libRadtran, which
%   defines the gamma distribution using the effective radius and the alpha
%   parameter, where alpha = 1/ve - 3


%   (5) inSitu_campaign - string containing the name of the field campaign

%   (6) inSitu_date - a datetime string containing the date of the sampled
%   data

%   (7) inSitu_time - double containing the UTC time of the samples half
%   way though the cloud (in altitude)


%   (8) compute_weighting_functions - a true or false flag that tells the
%   compute to compute weighting functions by creating N wc files where N
%   is equal to the number of cloud layers. The first file has all N
%   layers, and each successive file removes one cloud layer from the
%   bottom

%   (9) computer_name - the computer this code is running on

%   (10) index - this is the unique identifier that ensures files are not
%   written over one another. If one file is created, the index will be 1.
%   If many fiels are written in a for loop, each file will be tagged with
%   the number in the loop.


%   (11) wc_folder_path - Used to define where the water cloud wc files are
%   store




% OUTPUTS:
%   (1) .Dat file saved in the libRadTran folder:
%   /.../libRadtran-2.0.4/data/wc

% All look up tables were computed using the opensoure radiative transfer
% code, LibRadTran

% By Andrew John Buggee

%%

function [fileName] = write_wc_file_from_in_situ_2(re, lwc, z, alpha_width,...
    inSitu_campaign, inSitu_date, inSitu_time, compute_weighting_functions,...
    computer_name, index, wc_folder_path)

% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 11 inputs

if nargin~=10
    error([newline,'Should be 11 inputs: droplet effective radius, LWC, altitude,',...
        ' in-situ campaign name, date of data recording, time of data recording',...
        ', a flag telling the code to compute weighting functions, the computer name,',...
        'a unique file index,'...
        ' and the location where wc files will be saved.', newline])
end

% Check to make sure re is the same length as the altitude vector and the
% LWC vector

if length(re) ~= length(lwc) || length(re) ~= length(z) || length(re) ~= length(alpha_width)

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
    alpha_width = flipud(alpha_width);

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





