% Function to alter the above cloud column amount of the water vapor
% density profile

% INPUTS:

% (1) inputs - the input structure for creating the input files

% (2) aboveCloudTotal - This is the total column water vapor amount above
% cloud in kg/m^2

% (3) atm_folder_path - the location of libRadtrans atmosphere files and
% where to save the water vapor file

% By Andrew John Buggee

%%

function [filename_fullPath] = alter_aboveCloud_columnWaterVapor_profile(inputs, aboveCloudTotal, atm_folder_path)




% create the modified atm file name
fileName = [inputs.RT.atm_file(1:end-4), '_H2O_prof_', num2str(aboveCloudTotal), 'mm-aboveCloud.DAT'];



%% Read in the standard atmosphere profile

atm = read_libRadtran_atm_dat_profiles([atm_folder_path, inputs.RT.atm_file]);

% altitude values are listed in the 1st column (km)
z = atm(:,1);
% convert to kilometers
z = z*1e3; % (m)

% Water vapor densities are listed in the 7th column (cm^(-3))
waterVapor_column = atm(:, 7);
% convert water vapor densities to m^(-3)
waterVapor_column = waterVapor_column * 1e6;  % molecules/m^3


%% Solve for the scalar value that alters the above cloud column water vapor amount

% First, interpolate the profile so that the cloud top height and the
% sensor altitude are included
if ischar(inputs.RT.sensor_altitude) && strcmp(inputs.RT.sensor_altitude, 'toa')==true

    new_z = sort([z; inputs.RT.z_topBottom(1)*1e3], 'descend');  % m

    waterVapor_column_interp = interp1(z, waterVapor_column, new_z, "linear");    % moleules/m^3

    % set the index for the sensor location
    idx_sensor = 1;

elseif isdouble(inputs.RT.sensor_altitude)

    new_z = sort([z; inputs.RT.z_topBottom(1); inputs.RT.sensor_altitude], 'descending');

    waterVapor_column_interp = interp1(z, waterVapor_column, new_z, "linear");    % moleules/m^3

    % set the index for the sensor location
    idx_sensor = find(new_z==inputs.RT.sensor_altitude);

end
% Solve for the scalar value by integrating column water vapor from cloud
% top to sensor location

idx_cloudTop = find(new_z==(inputs.RT.z_topBottom(1)*1e3));

% integrate from cloud top to sensor location and convert the density
% profile from molecules/m^3 to kg/m^2
con = physical_constants;
aboveCloud_columnAmount = -(con.Mol_mass_h20_vap/con.N_A) *...
    trapz(new_z(idx_sensor:idx_cloudTop), waterVapor_column_interp(idx_sensor:idx_cloudTop));  % kg/m^2

% solve for the scalar constant
a = aboveCloudTotal / aboveCloud_columnAmount;


%% rescale the water vapor profile using the new scalar constant

waterVapor_column_2Write = waterVapor_column_interp;      % moleules/m^3

waterVapor_column_2Write(idx_sensor:idx_cloudTop) = a.*waterVapor_column_interp(idx_sensor:idx_cloudTop);  % moleules/m^3

% convert this back to molcules per cubic centimeter
waterVapor_column_2Write = waterVapor_column_2Write ./ 1e6;  % moleules/cm^3

% Convert the z vector to km
z_2write = new_z./1e3;  % (km)

%% Write the file

% Create the denisty file
fileID = fopen([atm_folder_path,fileName], 'w');

% fprintf writes lines in our text file from top to botom
% wc.DAT files are written with the higher altitudes at the top, and the
% surface at the bottom

% to write column vectors in a text file, we have to store them as row
% vectors
toWrite = [z_2write'; waterVapor_column_2Write'];

% Create the opening comment lines of the water vapor DAT file

fprintf(fileID, '%s %s %s \n','#','z', 'water-vapor-density');
fprintf(fileID, '%s %s %s \n','#','(km)','(molecules/cm^3)');

% Write in the data
fprintf(fileID,'%f %f \n', toWrite);
fclose(fileID);



%% define the output filename with full path

filename_fullPath = [atm_folder_path, fileName];

end
