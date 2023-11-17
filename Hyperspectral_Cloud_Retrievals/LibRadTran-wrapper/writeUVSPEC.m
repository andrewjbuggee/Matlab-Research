%% ----- WRITE UVSPEC INPUT FILES -----

% include a flag for verbose versus quiet, which will either create an
% error message text file or not




% By Andrew J. Buggee
%% Below are the instructions of how the file should be written

%   1) The very first line should be the rte_solver type. All values and
%   comments should be separated by only 1 space. This allows matlab to
%   simply read input settings.

%   2) 



function [outputArg1,outputArg2] = writeUVSPEC(inputArg1,inputArg2)
%WRITEUVSPEC Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;






%% ---- SETTING UP RTE SOLVER ----

% need to set the type of RTE solver:   
%       1) twostr (Two-Stream Approximation) - this is a 1D,
%       pseudo-spherical solver that calculates irradiance and actinic
%       flux. It is a very fast approximation, but it cannot solve for
%       radiances. 

%       2) disort - a 1D plane-parralel solver that calculates irradiance,
%       actinic flux and radiance. This is the default solver of UVSPEC,
%       and when in doubt, use disort. 


%% ---- SETTING UP SOLAR SOURCE SPECTRUM ----

% There are many options for a solar source. The defualt is simply 'solar'
% and, is using the default REPTRAN band parameterization, you can select
% wavelengths from 120 to 5000 nm. The options are:

%       1) source solar -- The default source which allows for spectrum
%       values between 120 and 5000 nm. 

%       2) source solar ../data/solar_flux/atlas_plus_modtran -- This high
%       resolution source only allows for a wavelengths between 200 and
%       799.95nm. This file does differ slightly from option 1, though it is
%       unclear why they differ at the time of this writing


%% --- Geometry Configuration ----

%   1) sza - Solar Zenith Angle - There can only be one per file. The
%   default value is 0, which implies that the sun is directly over had,
%   and along the line that intercepts with the earths center of mass. This
%   means the sun is on the horizon with an angle of 90. UVSPEC purports a
%   measured radiance of 0 over the solar spectrum at this angle. Angles
%   greater than 90 and up to 180 consittute a sun that is below the
%   horizon. Values between [0,180].  Units: degrees

%   4) phi0 - Solar Azimuth Angle - There can only be one per file. Units:
%   degrees. Values are between [0,360]

%   2) umu - Cosine(Zenith Viewing Angle) - Cosine of the viewing angle
%   where a viewing angle of 0 is straight down into the Earth, rather than
%   straight up. According to the manual, umu>0 is looking downward
%   (e.g. a satellite). So what this might mean is when umu>0 the default
%   sensor height is at TOA. umu<0 is looking upward. To make a
%   vector, just include spaces in between each value. If you set the
%   sensor altitude to be toa, then only values where umu>0 will give
%   non-zero radiance values, since umu>0 signals uvspec to look down.
%   Looking back into space will cause errors. If your sensor is at an
%   altitude of 10km, then a umu=0 is looking horizontally, umu>0 is looking
%   down towards the surface, and umu<0 is looking upwards towards the sky.
%   note: a umu of 0 is not allowed! This will lead to infinities because
%   1/umu when umu of 0 is infinity. And this is how we calculate the 
%   slant path. The values have to be increasing in order to be read. 
%   Units: degrees

%   3) phi - viewing aziumuth - this is the viewing azimuth of the sensor. This can
%   be input as a vector by leaving a blank space inbetween values. If a
%   vector is used, values must be increasing. Units: degrees. Values are
%   between [0,360], where 0 implies a device in the North looking
%   South, 90 implies a device in the East looking West, 180 implies a
%   device in the South looking North, and 270 implies a device in the West
%   looking East. 

%   5) zout - sensor altitude - the value is measured in km and must be
%   within the range of your atmospheric model defined by 'atmosphere
%   file'. For a satellite, set: zout TOA. This
%   will place the sensor at the top of the atmosphere. Then all umu values
%   should be greater than zero to specify looking downward. If zout is 0,
%   then umu should always be less than zero to specify that you are
%   looking upward. You can also set 'zout sur', to set the sensor at
%   the lowest atmospheric level. More than one value is allowed, but they
%   must be increasing in magnitude. Units: km

%   6) altitude - the level of the location above sea level, measured in
%   km. For example, if the input was 'altitude 1.25' then the lowest level
%   in the atmopshere model would be 1.25 km above sea level. The
%   atmosphere file is appropriately scaled to take this into account. Note
%   this is different from zout. If we set 'altitude 0.5' and 'zout 1' then
%   the surface is set at 0.5km above sea level and the sensor is 1 km
%   above the surface. Units: km

%% ---- Wavelength range and Slit Functions ----

%   1) slit_function_file - location of the slit function data file. This
%   for some reason has limitations. For example I cannot run a file with a
%   slit function over the range 400 - 800 nanometers, for whatever reason.
%   When I was using modtran we often used a slit function so the
%   transition between spectral bins was smooth. Atmospheric data taken
%   from measurements are located at discrete wavelengths. Slit functions
%   are used to interpolate across values and create more continuous
%   outputs. I should look into the different built in slit functions and
%   the limitations of each

%   2) wavelength - when only specifying a wavelength range using disort,
%   there are finite boundaries.


%% ---- Additional Output Quantities ----


% there is a host of additional outputs the user can opt for

%       1) output_quantity transmittance -- This outputs the spectral
%       transmittance of reflectance in place of all absolute quantities. 
%       Important note: the transmittance does not account for the 
%       sun-earth distance. It simply calculates I/I0. But
%       the transmittance DOES change for different SZA. I'm not sure why I
%       would want to use the transmissivity

%       2) output_quantity transmissivity -- The output calculates the
%       spectral transmittance instead of absolute quantities. But its
%       formula is different from option 1: I/(I0 * cos(x)) where
%       x is the solar zenith angle

%% ---- RUNNING CUSTOM FILES -----

% It appears there is a limit to the level of organizaiton LibRadTran will
% allow. Within the folder 'libRadtran-2.0.4', I'm only able to creat 1
% single additional folder. For example, I could make the directory
% '.../libRadtran-2.0.4/clouds' and here I can save several .INP file to
% run using uvspec. But if I added an additional directory, such as
% '.../libRadtran-2.0.4/clouds/opticallyThin', uvspec throws out an error
% saying it cannot find data files to create an atmosphere or wavelength
% grid. Apparently the logic in uvspec expects the directory 'data' to be
% only one level up from the current folder. This is stupid, but be
% careful!


end

