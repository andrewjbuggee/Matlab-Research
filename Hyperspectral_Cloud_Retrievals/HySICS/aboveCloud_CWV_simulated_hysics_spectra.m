% Determine above cloud column water vapor amount from simulated HySICS
% measurements


% INPUTS:

% (1) inputs - the input structure used to create the simulated
% measurements


% By Andrew John Buggee

%%

function [above_cloud_cwv] = aboveCloud_CWV_simulated_hysics_spectra(inputs)


% what computer is this function running on?
which_computer = whatComputer;

if strcmp(which_computer,'anbu8374')==true

    atm_folder_path = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/';

elseif strcmp(which_computer,'andrewbuggee')==true

    atm_folder_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
                       'LibRadTran/libRadtran-2.0.4/data/atmmod/'];
elseif strcmp(which_computer,'curc')==true

    atm_folder_path = '/projects/anbu8374/software/libRadtran-2.0.5/data/atmmod/';

end




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


%% Determine if the total column water vapor amount was altered

% define constants
con = physical_constants;

if inputs.RT.modify_total_columnWaterVapor==true

    % if true, then we must scale the profile accordingly
    c = inputs.RT.waterVapor_column / (-(con.Mol_mass_h20_vap/con.N_A) *trapz(z, waterVapor_column));

    % the new water vapor column is scaled 
    waterVapor_column = c.* waterVapor_column;

end

%% Compute the total column water vapor above cloud

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

above_cloud_cwv = -(con.Mol_mass_h20_vap/con.N_A) *...
    trapz(new_z(idx_sensor:idx_cloudTop), waterVapor_column_interp(idx_sensor:idx_cloudTop));  % kg/m^2





end
