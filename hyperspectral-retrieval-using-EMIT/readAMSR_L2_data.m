%% ----- READ IN AMSR DATA -----

% This function reads in AMSR-E and AMSR2 unified L2B ocean data as .h5 files
% Produces all necessary information in a structure

% Key fields extracted:
%   - Liquid Water Path (LiquidWaterPath)
%   - Total Precipitable Water (TotalPrecipitableWater)
%   - Wind Speed (WindSpeed)
%   - Sea Surface Temperature (ReynoldsSST)
%   - Geolocation (Latitude, Longitude)
%   - Quality flags and errors

% By Andrew J. Buggee
%%

function [AMSR] = readAMSR_L2_data(folderName)


% --- Add foldername to the matlab path ---
addpath(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName, '*.he5']);       % find all files that end in .h5

% check to see if we found any files!
if isempty(files)==true

    
    % check another variation of HDF5 files
    files = dir([folderName, '*.h5']);       % find all files that end in .h5

    if isempty(files)==true

    error([newline,'There are no .h5 files in the folder provided!', newline])

    end

end


for ii = 1:length(files)

    % Check if this is an AMSR file
    if contains(files(ii).name, 'AMSR') || contains(files(ii).name, 'amsr')

        fileName = [folderName, files(ii).name];
        
        %% ---- Read in HDF5 info structure -----
        info = h5info(fileName);

        %% ---- Read in Geolocation Fields ----
        % These are in the Geolocation_Fields group

        % Latitude and Longitude (npix x nscans)
        AMSR.geo.Latitude = double(h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Geolocation Fields/Latitude'));
        AMSR.geo.Longitude = double(h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Geolocation Fields/Longitude'));

        % Time for each scan (nscans)
        % Time is in seconds since 1993-01-01
        AMSR.Time = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Geolocation Fields/Time');

        % Convert time to UTC datetime
        % AMSR time starts on 1 Jan 1993
        epoch = datetime(1993,1,1,'TimeZone','UTC');
        AMSR.Time_UTC = epoch + seconds(AMSR.Time);
        
        % Read the human-readable time array (ntime x nscans)
        % This contains [year, month, day, hour, minute, second]
        AMSR.TimeHR = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/TimeHR');


        %% ---- Read in Data Fields ----
        % These are in the Data_Fields group

        % ------ Water Vapor ------
        % Total precipitable water (integrated water vapor) (npix x nscans)
        % units: mm
        AMSR.H2O.TotalPrecipitableWater = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/TotalPrecipitableWater');
        AMSR.H2O.ErrorTPW = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/ErrorTPW');

        % ------ Cloud Liquid Water ------
        % Liquid water path (cloud water only) (npix x nscans)
        % units: g/m^2
        AMSR.cloud.LiquidWaterPath = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/LiquidWaterPath');
        AMSR.cloud.ErrorLWP = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/ErrorLWP');

        % ------ Wind Speed ------
        % 10m wind speed (npix x nscans)
        % units: m/s
        AMSR.wind.WindSpeed = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/WindSpeed');
        AMSR.wind.ErrorWind = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/ErrorWind');

        % ------ Sea Surface Temperature ------
        % Reynolds SST (npix x nscans)
        % units: K
        AMSR.SST.ReynoldsSST = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/ReynoldsSST');

        % ------ Quality Control ------
        % Quality flag (npix x nscans)
        % 0: Highest quality retrieval
        % 1: Convergence reached
        % 2: No convergence, precipitation or land contamination possible
        % 3: TPW quality check failed
        % 4: Sun glint angle < 20 degrees (set to missing)
        % 5: Not run due to land or sea ice
        AMSR.QualityFlag = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/QualityFlag');

        % Chi-squared convergence metric (npix x nscans)
        AMSR.ChiSquared = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/ChiSquared');

        % % ------ Viewing Geometry ------
        % % Sun glint angle (npix x nscans)
        % AMSR.solar.GlintAngle = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/SunGlintAngle');
        % 
        % % ------ Surface Type ------
        % % Land percentage within FOV (npix x nscans)
        % AMSR.LandPercentage = h5read(fileName, '/HDFEOS/SWATHS/AMSR2_Level2_Ocean_Suite/Data Fields/LandPercentage');


        %% ---- Read in Global Attributes (Metadata) ----
        % Try to read global attributes for metadata
        try
            % Get the number of attributes
            num_attrs = length(info.Attributes);
            
            % Read all global attributes
            for jj = 1:num_attrs
                attr_name = info.Attributes(jj).Name;
                attr_value = info.Attributes(jj).Value;
                
                % Store in metadata structure
                % Replace spaces and special characters with underscores
                attr_name_clean = strrep(attr_name, ' ', '_');
                attr_name_clean = strrep(attr_name_clean, '-', '_');
                
                AMSR.metadata.(attr_name_clean) = attr_value;
            end
        catch
            warning('Could not read some global attributes');
        end


        %% ---- Apply Fill Value Masks ----
        % Replace fill values with NaN for easier data handling

        fillValue_float = -9999.0;
        fillValue_byte = -99;
        fillValue_short = -88;
        
        % Special marker values (according to comments in variable descriptions)
        landBadPixel = -998.0;  % land/bad pixel marker
        qualityIssue = -997.0;  % quality issue
        valueTooHigh = -996.0;  % value too high (TPW only)

        % Geolocation
        AMSR.geo.Latitude(AMSR.geo.Latitude == fillValue_float) = NaN;
        AMSR.geo.Longitude(AMSR.geo.Longitude == fillValue_float) = NaN;

        % Water vapor
        AMSR.H2O.TotalPrecipitableWater(AMSR.H2O.TotalPrecipitableWater == fillValue_float) = NaN;
        AMSR.H2O.TotalPrecipitableWater(AMSR.H2O.TotalPrecipitableWater == landBadPixel) = NaN;
        AMSR.H2O.TotalPrecipitableWater(AMSR.H2O.TotalPrecipitableWater == qualityIssue) = NaN;
        AMSR.H2O.TotalPrecipitableWater(AMSR.H2O.TotalPrecipitableWater == valueTooHigh) = NaN;
        AMSR.H2O.ErrorTPW(AMSR.H2O.ErrorTPW == fillValue_float) = NaN;

        % Cloud liquid water
        AMSR.cloud.LiquidWaterPath(AMSR.cloud.LiquidWaterPath == fillValue_float) = NaN;
        AMSR.cloud.LiquidWaterPath(AMSR.cloud.LiquidWaterPath == landBadPixel) = NaN;
        AMSR.cloud.LiquidWaterPath(AMSR.cloud.LiquidWaterPath == qualityIssue) = NaN;
        AMSR.cloud.ErrorLWP(AMSR.cloud.ErrorLWP == fillValue_float) = NaN;

        % Wind speed
        AMSR.wind.WindSpeed(AMSR.wind.WindSpeed == fillValue_float) = NaN;
        AMSR.wind.WindSpeed(AMSR.wind.WindSpeed == landBadPixel) = NaN;
        AMSR.wind.WindSpeed(AMSR.wind.WindSpeed == qualityIssue) = NaN;
        AMSR.wind.ErrorWind(AMSR.wind.ErrorWind == fillValue_float) = NaN;

        % Sea surface temperature
        AMSR.SST.ReynoldsSST(AMSR.SST.ReynoldsSST == fillValue_float) = NaN;

        % Chi-squared
        AMSR.ChiSquared(AMSR.ChiSquared == fillValue_float) = NaN;
        AMSR.ChiSquared(AMSR.ChiSquared == landBadPixel) = NaN;
        AMSR.ChiSquared(AMSR.ChiSquared == qualityIssue) = NaN;

        % Sun glint angle
        % AMSR.solar.GlintAngle(AMSR.solar.GlintAngle == fillValue_short) = NaN;

        % Quality flag (keep as-is, it's categorical)
        % Land percentage (keep as-is, 0-100%)


        %% ---- Display Summary ----
        % fprintf('AMSR data successfully loaded!\n');
        % fprintf('File: %s\n', files(ii).name);
        % fprintf('Dimensions: %d pixels x %d scans\n', size(AMSR.geo.Latitude,1), size(AMSR.geo.Latitude,2));
        % fprintf('Latitude range: [%.2f, %.2f]\n', min(AMSR.geo.Latitude(:), [], 'omitnan'), max(AMSR.geo.Latitude(:), [], 'omitnan'));
        % fprintf('Longitude range: [%.2f, %.2f]\n', min(AMSR.geo.Longitude(:), [], 'omitnan'), max(AMSR.geo.Longitude(:), [], 'omitnan'));
        % 
        % % Display data ranges
        % fprintf('\nData ranges:\n');
        % fprintf('  TPW: [%.2f, %.2f] mm\n', min(AMSR.H2O.TotalPrecipitableWater(:), [], 'omitnan'), ...
        %     max(AMSR.H2O.TotalPrecipitableWater(:), [], 'omitnan'));
        % fprintf('  LWP: [%.2f, %.2f] g/m^2\n', min(AMSR.cloud.LiquidWaterPath(:), [], 'omitnan'), ...
        %     max(AMSR.cloud.LiquidWaterPath(:), [], 'omitnan'));
        % fprintf('  Wind Speed: [%.2f, %.2f] m/s\n', min(AMSR.wind.WindSpeed(:), [], 'omitnan'), ...
        %     max(AMSR.wind.WindSpeed(:), [], 'omitnan'));
        % fprintf('  SST: [%.2f, %.2f] K\n', min(AMSR.SST.ReynoldsSST(:), [], 'omitnan'), ...
        %     max(AMSR.SST.ReynoldsSST(:), [], 'omitnan'));


    end

end


end