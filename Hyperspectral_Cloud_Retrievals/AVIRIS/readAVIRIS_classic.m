%% ----- READ IN AVIRIS DATA -----


%   2) dataType - 1= Byte: 8-bit unsigned integer (unit8); 2= Integer: 16-bit signed integer (int16);
%   3= Long: 32-bit signed integer (int32); 4=Floating-Point: 32-bit
%   single-precision (single); 5= Double-Precision: 64-bit double precision
%   floating point (double); 6= Complex: Read-Imaginary pair of signle-precision
%   floating-point (complex)

%   3) headerOffset - The number of bytes that make up the header of a
%   given file. Often the header offset is 0, and the header information is
%   simply stored in a seperate file altogether.

%   4) interLeave - This tells matlab how the band data is interleaved, or
%   stored, in the data cube. The options are 1) 'bsq' - or band
%   sequential, which is the simplest format. 2) 'bip' - or band
%   interleaved by pixel. 3) 'bil' - or band interleaved by line.

%   5) byteOrder - Order of bytes in integer, long integer, 64-bit integer,
%   unsigned 64-bit integer, floating point, double precision and complex.
%   0 = 'native' (Host); 1 = 'ieee-be' (Netowrk)


% The ouput being extracted is calibrated radiance, which has units of
% micro-watts per centimeter squared per nanometer per steradian



% By Andrew J. Buggee
%% --- DEFINE FOLDER AND FILE NAMES ---

function [aviris] = readAVIRIS_classic(folderName, justPosition)

% folderName = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/AVIRIS/AVIRIS-Classic-data/8_20_2018/orthorectified_radiance_data/'];

% fileName = 'f180820t01p00r09rdn_e';

%% --- READ THE HEADER FILE ---

% find radiance files
radiance_hdr_files = dir([folderName,'*ort_img.hdr']);
radiance_files = dir([folderName,'*_ort_img']);

% The gain file is needed to convert radiance Int16 to double
gain_file = dir([folderName,'*.gain']);

% find sensor position data
position_files = dir([folderName, '*ort_igm']);
position_files_hdr = dir([folderName, '*ort_igm.hdr']);

% find observation geometry data that matches the corrected orthorectified
% data
observation_geometry_hdr = dir([folderName, '*b_obs_ort.hdr']);
observation_geometry = dir([folderName, '*b_obs_ort']);



if justPosition==true
    
    
    % ---------------------------------------------------------
    % ------------ Extract AVIRIS Position Data ---------------
    % ---------------------------------------------------------
    
    % The AVIRIS position data is stored in the Band-Interleaved-by-Pixel
    % format
    % These are the parameters relating to the geometry of the observation
    
    %   (1) longitude (degrees) - from -180 to 180. values less than 0 are
    %   measured in degrees west of Greenwich England.
    %   (2) latitude (degrees) - from -90 to 90. Values greater than 0
    %   represent the northern hemisphere
    %   (3) ground elevation at each pixel (meters)
    
    info = enviinfo([folderName,position_files_hdr.name]);
    
    
    dataDim = [info.Height,info.Width,info.Bands]; % lines, samples, bands. taken from obs hdr file
    dataType = info.DataType; % taken from obs hdr file
    headerOffset = info.HeaderOffset; % header offset. taken from obs hdr file
    interLeave = info.Interleave; % taken from obs hdr file
    byteOrder = info.ByteOrder; % taken from obs hdr file
    
    
    
    aviris.position = multibandread([folderName,position_files.name],dataDim,dataType,headerOffset,interLeave,byteOrder);
    
else
    
    
    % ---------------------------------------------------------
    % --------------- Extract radiance data -------------------
    % ---------------------------------------------------------
    
    % The radiance data is stored as Band Interleaved by Pixel (BIP)
    
    info = enviinfo([folderName,radiance_hdr_files.name]);
    gain = importdata([folderName, gain_file.name]); % import gain values that convert int16 to radiance
    
    
    dataDim = [info.Height,info.Width,info.Bands]; % lines, samples, bands. taken from _ort_img.hdr file
    dataType = info.DataType; % Taken from _obs.hdr file
    headerOffset = info.HeaderOffset; % header offset. Taken from _obs.hdr file
    interLeave = info.Interleave; % Taken from _obs.hdr file
    byteOrder = info.ByteOrder; % Taken from _obs.hdr file
    
    
    
    aviris.radiance = multibandread([folderName,radiance_files.name],dataDim,dataType,headerOffset,interLeave,byteOrder); % mu-W/cm^2/nm/sr
    
    % convert radiance from int16 to micro-watts/cm^2/nm/sr. Accroding to the
    % readme file: "When each spectrum is divided by the factors in this file
    % the 16-bit integers are converted to radiance in units of
    % (micro-watts/cm^2/nm/sr)"
    
    gain_vals = reshape(gain(:,1),1,1,[]);
    aviris.radiance = aviris.radiance./repmat(gain_vals,size(aviris.radiance,1),size(aviris.radiance,2),1);  % micro-watts/cm^2/nm/sr
    
    
    
    % ---------------------------------------------------------
    % ------------ Extract Observation Geometry ---------------
    % ---------------------------------------------------------
    
    % The observation geometry data is stored as Band Interleaved by Pixel (BIP)
    % These are the parameters relating to the geometry of the observation
    % The 10 values provided are:
    %
    %   (1) path length (sensor to ground) in meters
    %   (2) to-sensor-azimuth (0 to 360 degrees clockwise from North)
    %   (3) to-sensor-zenith (0 to 90 degrees zenith)
    %   (4) to-sun-azimuth (0 to 360 degrees clockwise from North)
    %   (5) to-sun-zenith (0 to 90 degrees zenith)
    %   (6) solar phase (degrees between to-sensor and to-sun vectors in
    %   principal plane)
    %   (7) slope (local surface slope as derived from DEM in degrees)
    %   (8) aspect (local surface aspect 0 to 360 degrees clockwise from N)
    %   (9) cosine i (apparent local illumination factor based on DEM slop and
    %   aspect and to-sun vector. Ranges from -1 to 1)
    %   (10) UTC time (decimal hours for mid-line pixels
    
    info = enviinfo([folderName,observation_geometry_hdr.name]);
    
    
    dataDim = [info.Height,info.Width,info.Bands]; % lines, samples, bands. taken from obs hdr file
    dataType = info.DataType; % taken from obs hdr file
    headerOffset = info.HeaderOffset; % header offset. taken from obs hdr file
    interLeave = info.Interleave; % taken from obs hdr file
    byteOrder = info.ByteOrder; % taken from obs hdr file
    
    
    
    aviris.obsGeometry = multibandread([folderName,observation_geometry.name],dataDim,dataType,headerOffset,interLeave,byteOrder);
    
    
    
    
    
    % ---------------------------------------------------------
    % ------------ Extract AVIRIS Position Data ---------------
    % ---------------------------------------------------------
    
    % The AVIRIS position data is stored in the Band-Interleaved-by-Pixel
    % format
    % These are the parameters relating to the geometry of the observation
    
    %   (1) longitude (degrees
    %   (2) latitude (degrees
    %   (3) ground elevation at each pixel (meters)
    
    info = enviinfo([folderName,position_files_hdr.name]);
    
    
    dataDim = [info.Height,info.Width,info.Bands]; % lines, samples, bands. taken from obs hdr file
    dataType = info.DataType; % taken from obs hdr file
    headerOffset = info.HeaderOffset; % header offset. taken from obs hdr file
    interLeave = info.Interleave; % taken from obs hdr file
    byteOrder = info.ByteOrder; % taken from obs hdr file
    
    
    
    aviris.position = multibandread([folderName,position_files.name],dataDim,dataType,headerOffset,interLeave,byteOrder);
    
    
end




end

