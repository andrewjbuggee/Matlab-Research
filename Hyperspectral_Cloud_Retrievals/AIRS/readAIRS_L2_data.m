%% ----- READ IN AIRS DATA -----

% This function reads in L2 AIRS data as .hdf files
% Produces all necessary information in a structure

% Key fields extracted:
%   - Atmospheric temperature profiles (TAirStd)
%   - Pressure levels (pressStd)
%   - Water vapor profiles (H2OMMRStd)
%   - Methane profiles (CH4VMRLevStd)
%   - Carbon monoxide profiles (COVMRLevStd)
%   - Geolocation (Latitude, Longitude)
%   - Quality control flags

% By Andrew J. Buggee
%%

function [airs] = readAIRS_L2_data(folderName)


% --- Add foldername to the matlab path ---
addpath(folderName)

% -----------------------------------------------------
% ----- Check to see if the folder name is valid ------
% -----------------------------------------------------


files = dir([folderName, '*.hdf']);       % find all files that end in .hdf

% check to see if we found any files!
if isempty(files)==true
    error([newline,'There are no files in the folder provided!', newline])
end


for ii = 1:length(files)



    if strcmp(files(ii).name(1:4), 'AIRS')==true

        fileName = [folderName, files(ii).name];
        %% ---- Read in HDF info structure -----
        info = hdfinfo(fileName);

        %% ---- Read in Geolocation Fields ----
        % These are in the Geolocation_Fields group

        % Latitude and Longitude (GeoTrack x GeoXTrack)
        airs.geo.Latitude = hdfread(fileName, 'Latitude');
        airs.geo.Longitude = hdfread(fileName, 'Longitude');

        % Time for each footprint (seconds since 1993-01-01, TAI calendar)
        airs.Time = hdfread(fileName, 'Time');

        % Convert time to UTC datetime
        % AIRS time starts on 1 Jan 1993 (TAI)
        epoch = datetime(1993,1,1,'TimeZone','UTC');
        airs.Time_UTC = epoch + seconds(airs.Time);


        %% ---- Read in Data Fields ----
        % These are in the Data_Fields group

        % ------ Pressure Levels ------
        % Standard pressure levels for temperature and other retrievals
        airs.pressStd = hdfread(fileName, 'pressStd');  % (StdPressureLev = 28)   % hPa

        % Pressure levels for water vapor retrievals
        airs.pressH2O = hdfread(fileName, 'pressH2O');  % (H2OPressureLev = 15)  % hpa


        % ------ Temperature Profiles ------
        % Air temperature at standard pressure levels (GeoTrack x GeoXTrack x StdPressureLev)
        airs.temp.prof_std = hdfread(fileName, 'TAirStd');         % K
        airs.temp.prof_std_QC = hdfread(fileName, 'TAirStd_QC');
        airs.temp.prof_std_err = hdfread(fileName, 'TAirStdErr');

        % Surface air temperature
        airs.temp.SurfAir = hdfread(fileName, 'TSurfAir');       % K
        airs.temp.SurfAir_QC = hdfread(fileName, 'TSurfAir_QC');

        % Surface pressure
        airs.pressure.SurfStd = hdfread(fileName, 'PSurfStd');
        airs.pressure.SurfStd_QC = hdfread(fileName, 'PSurfStd_QC');


        % ------ Water Vapor Profiles ------
        % Water vapor mass mixing ratio at H2O pressure layers (GeoTrack x GeoXTrack x H2OPressureLay)
        airs.H2O.MMR_Std = hdfread(fileName, 'H2OMMRStd');
        airs.H2O.MMR_Std_QC = hdfread(fileName, 'H2OMMRStd_QC');
        airs.H2O.MMR_StdErr = hdfread(fileName, 'H2OMMRStdErr');

        % Water vapor at pressure levels (not layers)
        airs.H2O.MMR_LevStd = hdfread(fileName, 'H2OMMRLevStd');        % Layer water Vapor Mass Mixing Ratio (gm/kg dry air)  
        airs.H2O.MMR_LevStd_QC = hdfread(fileName, 'H2OMMRLevStd_QC');

        % Total column water vapor
        airs.H2O.totCol_Std = hdfread(fileName, 'totH2OStd');      % Total precipitable water vapor (kg/m2)
        airs.H2O.totCol_Std_QC = hdfread(fileName, 'totH2OStd_QC');

        % Relative humidity
        airs.H2O.RelHum = hdfread(fileName, 'RelHum');           % % RH
        airs.H2O.RelHum_QC = hdfread(fileName, 'RelHum_QC');


        % ------ Methane Profiles ------
        % Methane volume mixing ratio at standard pressure levels (GeoTrack x GeoXTrack x StdPressureLev)
        airs.CH4.VMR_LevStd = hdfread(fileName, 'CH4VMRLevStd');
        airs.CH4.VMR_LevStd_QC = hdfread(fileName, 'CH4VMRLevStd_QC');
        airs.CH4.VMR_LevStdErr = hdfread(fileName, 'CH4VMRLevStdErr');

        % Total column methane
        airs.CH4.total_column = hdfread(fileName, 'CH4_total_column');
        airs.CH4.total_column_QC = hdfread(fileName, 'CH4_total_column_QC');

        % Degrees of freedom for methane retrieval
        airs.CH4.dof = hdfread(fileName, 'CH4_dof');


        % ------ Carbon Monoxide Profiles ------
        % CO volume mixing ratio at standard pressure levels (GeoTrack x GeoXTrack x StdPressureLev)
        airs.CO.VMR_LevStd = hdfread(fileName, 'COVMRLevStd');
        airs.CO.VMR_LevStd_QC = hdfread(fileName, 'COVMRLevStd_QC');
        airs.CO.VMR_LevStdErr = hdfread(fileName, 'COVMRLevStdErr');

        % Total column CO
        airs.CO.total_column = hdfread(fileName, 'CO_total_column');
        airs.CO.total_column_QC = hdfread(fileName, 'CO_total_column_QC');

        % Degrees of freedom for CO retrieval
        % AIRS.CO_dof = hdfread(fileName, 'CO_dof');


        % ------ Ozone Profiles (bonus, might be useful) ------
        % Ozone volume mixing ratio
        airs.O3.VMR_Std = hdfread(fileName, 'O3VMRStd');
        airs.O3.VMR_Std_QC = hdfread(fileName, 'O3VMRStd_QC');

        % Total column ozone
        airs.O3.totCol_Std = hdfread(fileName, 'totO3Std');
        airs.O3.totCol_Std_QC = hdfread(fileName, 'totO3Std_QC');


        % ------ Cloud Properties ------
        % Total cloud fraction
        airs.cloud.FracTot = hdfread(fileName, 'CldFrcTot');
        airs.cloud.FracTot_QC = hdfread(fileName, 'CldFrcTot_QC');

        % Cloud top pressure
        airs.cloud.presTop_ = hdfread(fileName, 'PCldTop');
        airs.cloud.presTop_QC = hdfread(fileName, 'PCldTop_QC');


        % ------ Viewing Geometry ------
        % Satellite zenith angle
        airs.sensor.zen = hdfread(fileName, 'satzen');

        % Solar zenith angle
        airs.solar.zen = hdfread(fileName, 'solzen');

        % Surface topography
        % AIRS.topog = hdfread(fileName, 'topog');

        % Land fraction (0 = ocean, 1 = land)
        airs.landFrac = hdfread(fileName, 'landFrac');


        %% ---- Read in Swath Attributes (Metadata) ----
        % These contain useful granule-level information

        % Try to read attributes - structure varies by HDF implementation
        try
            airs.metadata.processing_level = hdfread(fileName, 'processing_level');
            airs.metadata.DayNightFlag = hdfread(fileName, 'DayNightFlag');
            airs.metadata.start_year = hdfread(fileName, 'start_year');
            airs.metadata.start_month = hdfread(fileName, 'start_month');
            airs.metadata.start_day = hdfread(fileName, 'start_day');
            airs.metadata.start_hour = hdfread(fileName, 'start_hour');
            airs.metadata.start_minute = hdfread(fileName, 'start_minute');
            airs.metadata.granule_number = hdfread(fileName, 'granule_number');
        catch
            warning('Could not read some metadata attributes from Swath_Attributes group');
        end


        %% ---- Apply Fill Value Masks ----
        % Replace fill values with NaN for easier data handling

        fillValue_float = -9999.0;
        fillValue_short = -9999;

        % Temperature
        airs.temp.prof_std(airs.temp.prof_std == fillValue_float) = NaN;
        airs.temp.SurfAir(airs.temp.SurfAir == fillValue_float) = NaN;

        % Pressure
        airs.pressure.SurfStd(airs.pressure.SurfStd == fillValue_float) = NaN;

        % Water vapor
        airs.H2O.MMR_Std(airs.H2O.MMR_Std == fillValue_float) = NaN;
        airs.H2O.MMR_LevStd(airs.H2O.MMR_LevStd == fillValue_float) = NaN;
        airs.H2O.totCol_Std(airs.H2O.totCol_Std == fillValue_float) = NaN;
        airs.H2O.RelHum(airs.H2O.RelHum == fillValue_float) = NaN;

        % Methane
        airs.CH4.VMR_LevStd(airs.CH4.VMR_LevStd == fillValue_float) = NaN;
        airs.CH4.total_column(airs.CH4.total_column == fillValue_float) = NaN;

        % Carbon monoxide
        airs.CO.VMR_LevStd(airs.CO.VMR_LevStd == fillValue_float) = NaN;
        airs.CO.total_column(airs.CO.total_column == fillValue_float) = NaN;

        % Ozone
        airs.O3.VMR_Std(airs.O3.VMR_Std == fillValue_float) = NaN;
        airs.O3.totCol_Std(airs.O3.totCol_Std == fillValue_float) = NaN;

        % Clouds
        airs.cloud.FracTot(airs.cloud.FracTot == fillValue_float) = NaN;
        airs.cloud.presTop_(airs.cloud.presTop_ == fillValue_float) = NaN;

        % Geometry
        airs.sensor.zen(airs.sensor.zen == fillValue_float) = NaN;
        airs.solar.zen(airs.solar.zen == fillValue_float) = NaN;

        % AIRS.topog(AIRS.topog == fillValue_float) = NaN;
        airs.landFrac(airs.landFrac == fillValue_float) = NaN;


        %% ---- Display Summary ----
        % fprintf('AIRS data successfully loaded!\n');
        % fprintf('Dimensions: %d GeoTrack x %d GeoXTrack\n', size(airs.geo.Latitude,1), size(airs.geo.Latitude,2));
        % fprintf('Temperature levels: %d\n', length(airs.pressStd{1}));
        % fprintf('Water vapor levels: %d\n', length(airs.pressH2O{1}));
        % fprintf('Latitude range: [%.2f, %.2f]\n', min(airs.geo.Latitude(:)), max(airs.geo.Latitude(:)));
        % fprintf('Longitude range: [%.2f, %.2f]\n', min(airs.geo.Longitude(:)), max(airs.geo.Longitude(:)));


    end

end


end