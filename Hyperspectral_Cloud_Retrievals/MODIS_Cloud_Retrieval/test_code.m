%% Read in MODIS L1B Swath Meta Data

L1B_metadata = hdfread(fileName, 'Level 1B Swath Metadata');

raw_time = L1B_metadata{5};     % TAI time (number of sec. since 1/1/93)

% convert this to UTC time
% MODIS EV Sector time starts on 1 Jan 1993
epoch = datetime(1993,1,1,'TimeZone','UTCLeapSeconds');

% now lets add secods according to the EV sector start time measurements
modis_pixel_time = epoch + seconds(raw_time);


%%

