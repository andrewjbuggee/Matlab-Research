%% Scratch work to figure out why my calcualtions of VOCALS Rex data differ from Painemal and Zudema

% By Andrew John Buggee

%%  Grab VOCALS-REx Data

clear variables

scriptPlotting_wht;

% define the VOCALS-REx data folder
foldername = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/';

% define filenames
% ----- November 9th data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ----- November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% load vocals rex data
vr = readVocalsRex([foldername, filename]);

%% Grab MODIS data

% ----- November 9th at decimal time 0.611 (14:40) -----
%modis_folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/';


% ----- November 11th at decimal time 0.604 (14:30) -----
%modis_folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/';


% ----- November 11th at decimal time 0.784 (18:50) -----
modis_folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/';


% ----- November 13th at decimal time 0.663 (15:55) -----
%modis_folderName = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_13/';



[modis,L1B_500m_fileName] = retrieveMODIS_data(modis_folderName);

modis_dataTime = (modis.time(1)*3600 + modis.time(2)*60)/86400;


%% Painemal and Zudema use decimal time

% Want to turn on diagnostic tools and look at some plots of MODIS
% estimates around the mimima value?

diagnostic_tools = false;

% =========================================================================
% *************************************************************************
% -------------------------------------------------------------------------

% YOU MUST MANUALLY ENTER THE TIME OF THE CLOUD PROFILE YOUR INTERESTED IN!


% This is used to identify when along a flights path they assumed a cloud
% and made their calculations

% Lets use their decimal indicator to find where in the data to compare

pz_time = 0.7816;


% -------------------------------------------------------------------------
% *************************************************************************
% =========================================================================



% Some data sets need to be aligned with a MODIS data set by using an index
% value a time time steps away from the PZ_time

if pz_time==0.6069 && strcmp(modis_folderName,'/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1430/')==true

    index_delay = 4;

elseif pz_time==0.7816 && strcmp(modis_folderName,'/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_11_1850/')==true

    index_delay = -8;


elseif pz_time==0.6120 && strcmp(modis_folderName,'/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/MODIS_Cloud_Retrieval/MODIS_data/2008_11_09/')==true

    index_delay = 0;

else

    index_delay = 0;

end


% convert the start time to decimal time
sec_per_hour = 3600;                                                                    % sec
sec_per_min = 60;                                                                       % sec
sec_per_day = 86400;                                                                    % sec

startTime_dec = (vr.startTime(1)*3600 + vr.startTime(2)*60)/86400;                      % fraction of a day 

pz_secSinceStart = floor((pz_time - startTime_dec)*sec_per_day);                               % seconds since flight started up to the time indicated by painemal and zudema

modis_secSinceStart = floor((modis_dataTime - startTime_dec)*sec_per_day);


% Define the number of discrete points in time you wish to look at. The
% center point will be at the start time calculated above
windowLength = 100;

% Using the time defined above, find the location in space closest to
% the airplanes location

% Using lat and long with can minimize the euclidean norm
% ***** There are some special cases *****

% The values are closer on Novemeber 11th 2008 if we compare one time step
% ahead. The re values are spatially similar, but the optical thickness
% change quite a bit!
dist_btwn_PZ_startTime_and_MODIS = sqrt((modis.geo.lat - vr.latitude(pz_secSinceStart+index_delay)).^2 + (modis.geo.long - vr.longitude(pz_secSinceStart+index_delay)).^2); 
[minDist, index_minDist] = min(dist_btwn_PZ_startTime_and_MODIS, [], 'all');


% Find the row and column locations of the minima
[row,col] = ind2sub(size(modis.cloud.effRadius17), index_minDist);

if diagnostic_tools == true
    
    
    x = modis.cloud.effRadius17(row-2:row+2, col-2:col+2);
    y = modis.cloud.optThickness17(row-2:row+2, col-2:col+2);
    c = dist_btwn_PZ_startTime_and_MODIS(row-2:row+2, col-2:col+2);

    figure; subplot(1,2,1)
    s = scatter(x(:),y(:),80, c(:), 'filled');
    grid on; grid minor;
    colorbar
    xlabel('Effective Radius (\mum)')
    ylabel('Optical Thickness')
    title('Bands 1 & 7 - Color = min(dist)')
    
    x = modis.geo.long(row-2:row+2, col-2:col+2);
    y = modis.geo.lat(row-2:row+2, col-2:col+2);



    subplot(1,2,2); s = scatter(x(:),y(:),80, c(:), 'filled');
    grid on; grid minor;
    colorbar
    xlabel('Longitude')
    ylabel('Latitude')
    title('Bands 1 & 7 - Color = min(dist)')
    hold on;

    set(gcf, 'Position', [0 0 1300 400])
    

end

% grab the value for droplet size at the location closest to the planes
% locaiton at the time specified above
re_16_modis = modis.cloud.effRadius16(index_minDist);
re_17_modis = modis.cloud.effRadius17(index_minDist);

% Grab the optical depth values

tau_16_modis = modis.cloud.optThickness16(index_minDist);
tau_17_modis = modis.cloud.optThickness17(index_minDist);

%% Make a plot of the number concentration and effective radius aat the time defined above

% ********* USING SUMMATION TO COMPUTE EFFECTIVE RADIUS *********

% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];

figure; subplot(1,2,1)
plot(vr.time, vr.re,'k'); xlabel('Time (sec)'); ylabel('r_{e} (\mum)')
grid on; grid minor; hold on; 

% Plot the modis droplet estimate as a constant vertical line
yl0 = yline(re_17_modis,'--',['MODIS $$r_{e} = $$',num2str(re_17_modis), '$$\mu m$$'], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
yl0.LabelVerticalAlignment = 'top';

% now plot the altitude of the plane on the right hand side
yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m)')
% zoom in on a particular window

xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);



% plot a vertical line marking the bottom of the cloud
% This occurs when the LWC is atleast 0.03 c/m^3
lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

% Plot a vertical line at the moment where LWC first exceeds the limit set
% above
LWC_window_index = pz_secSinceStart - windowLength/2 : pz_secSinceStart + windowLength/2;
LWC_window = vr.lwc(LWC_window_index);
index_minVal = find(LWC_window >= lwc_lim);

% Plot a vertical line!
xl = xline(vr.time(LWC_window_index(1) + index_minVal(1) -1),'--','$$LWC = 0.03\; g/m^{3}$$', 'Fontsize',18, 'Interpreter','latex','LineWidth',2);
xl.LabelVerticalAlignment = 'top';


subplot(1,2,2)
plot(vr.time, vr.total_Nc,'k'); xlabel('Time (sec)'); ylabel('N_{c} (cm^{-3})')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m)')
% zoom in on a particular window
xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

% Plot a vertical line!
xl = xline(vr.time(LWC_window_index(1) + index_minVal(1) -1),'--','$$LWC = 0.03\; g/m^{3}$$', 'Fontsize',18, 'Interpreter','latex','LineWidth',2);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = "left";

set(gcf, 'Position', [0 0 1300 400])


%% Plot droplet size vs liquid water content
% Using Painemal and Zudema definition for cloud top and bottom

lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

% I want a subplot with the number concentration and liquid water content,
% and the effective radius with liquid water content

figure; subplot(1,2,1)
plot(vr.time, vr.re); xlabel('Time (sec)'); ylabel('r_{e} (\mum)')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.lwc); ylabel('LWC (g/m^{3})')
% zoom in on a particular window

xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

% Plot a vertical line at the moment where LWC first exceeds the limit set
% above
LWC_window_index = pz_secSinceStart - windowLength/2 : pz_secSinceStart + windowLength/2;
LWC_window = vr.lwc(LWC_window_index);
index_minVal = find(LWC_window >= lwc_lim);

% Plot a vertical line!
xl = xline(vr.time(LWC_window_index(1) + index_minVal(1) -1),'--','$$LWC = 0.03\; g/m^{3}$$', 'Fontsize',18, 'Interpreter','latex','LineWidth',2);
xl.LabelVerticalAlignment = 'top';


subplot(1,2,2)
plot(vr.time, vr.total_Nc); xlabel('Time (sec)'); ylabel('N_{c} (cm^{-3})')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.lwc); ylabel('LWC (g/m^{3})')
% zoom in on a particular window
xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

set(gcf, 'Position', [0 0 1300 400])



%% ---- Make a plot of the effective radius within the boundary limits of the cloud ----


% ********* USING SUMMATION TO COMPUTE EFFECTIVE RADIUS *********
% ********* Compute optical depth using cloud boundaries defined by... *********

% Cloud bottom = location where LWC = 0.03 g/m^3
% ****** There are two ways we could define LWC *****
% Cloud top is tricky. The vocals rex data shows droplets being recorded
% for several data points past the peak LWC value. But MODIS estimates for
% optical depth better align with the cloud top being defined as the peak
% LWC value. After this maximum, LWC tends to fall off dramatically

% Cloud top = location where re goes to 0 after LWC> 0.03 g/m^3
% Cloud top = maximum valude of LWC after LWC> 0.03 g/m^3

% First find the index cloud bottom using the definition above

lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

LWC_window_index = pz_secSinceStart - windowLength/2 : pz_secSinceStart + windowLength/2;
LWC_window = vr.lwc(LWC_window_index);
index_minVal = find(LWC_window >= lwc_lim);
index_cloudBottom = LWC_window_index(1) + index_minVal(1) -1;



% Now lets find the index for cloud top using the definition above
% find the first value of re=0 after the above index

window_data = vr.re(index_cloudBottom : index_cloudBottom + 100);
index_nan = find(isnan(window_data));

index_cloudTop = index_cloudBottom + index_nan(1) - 2;

% Lets also find the index where LWC is a maxium after the index found for
% cloud bottom

window_data = vr.lwc(index_cloudBottom : index_cloudBottom + 100);
[~,index_maxLWC] = max(window_data);

index_cloudTop2 = index_cloudBottom + index_maxLWC - 1;


% ------------------------------------------------------------------
% ------------------ Compute optical depth -------------------------
% ------------------------------------------------------------------
% optical depth is defined to be 0 at cloud top and increasing towards
% cloud bottom
vector_length = length((index_cloudBottom:index_cloudTop));
tau = zeros(1,vector_length-1);
for ii = 1:length(vr.altitude(index_cloudBottom:index_cloudTop))-1
    
    % we have to convert Nc and re to have the same units as the alitude,
    % which is in meters
    re_meters = vr.re(index_cloudBottom+vector_length-1-ii:index_cloudTop)./1e6;                                   % meters
    total_Nc_meters = vr.total_Nc(index_cloudBottom+vector_length-1-ii:index_cloudTop).*1e6;                           % #/m^3
    altitude = vr.altitude(index_cloudBottom+vector_length-1-ii:index_cloudTop);
    % we need to flip these vectors so we start integratin at the cloud
    % top!
    tau(ii) = -2*pi* trapz(fliplr(altitude), fliplr(re_meters.^2 .* total_Nc_meters));

end

% add a zero at the begining!
tau = [0,tau];


% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];
nice_orange = [0.8500, 0.3250, 0.0980];

figure;
plot(fliplr(vr.re(index_cloudBottom:index_cloudTop)),tau,'k'); 
set(gca,'YDir','reverse')
ylabel('$$\tau$$','interpreter','latex'); xlabel('$$r_{e}$$ $$(\mu m)$$','Interpreter','latex')
title('Cloud Top = $$(r_e = 0)$$', 'Interpreter','latex')
grid on; grid minor; hold on; 

% Plot the modis droplet estimate as a constant vertical line
xl0 = xline(re_16_modis,'--',['MODIS $$r_{1.6} = $$',num2str(re_16_modis), '$$\mu m$$'], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(tau_16_modis,'--',['MODIS $$\tau_{1.6} = $$',num2str(tau_16_modis)], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';

% Let's compute the mean number concentration within this cloud and print
% it on our plot
mean_Nc = mean(vr.total_Nc(index_cloudBottom:index_cloudTop));

dim = [.2 .5 .3 .3];
str = ['$$< N_c > = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18,'FontWeight','bold');



% ----- Using the Second definition for cloud top ------

% ------------------------------------------------------------------
% ------------------ Compute optical depth -------------------------
% ------------------------------------------------------------------
% optical depth is defined to be 0 at cloud top and increasing towards
% cloud bottom
vector_length = length((index_cloudBottom:index_cloudTop2));
tau = zeros(1,vector_length-1);
for ii = 1:length(vr.altitude(index_cloudBottom:index_cloudTop2))-1
    
    % we have to convert Nc and re to have the same units as the alitude,
    % which is in meters
    re_meters = vr.re(index_cloudBottom+vector_length-1-ii:index_cloudTop2)./1e6;                                   % meters
    total_Nc_meters = vr.total_Nc(index_cloudBottom+vector_length-1-ii:index_cloudTop2).*1e6;                           % #/m^3
    %altitude = vr.altitude(index_cloudBottom+vector_length-1-ii:index_cloudTop2);
    altitude = vr.altitude(index_cloudTop2) -  vr.altitude(index_cloudBottom+vector_length-1-ii:index_cloudTop2);
    % we need to flip these vectors so we start integratin at the cloud
    % top!
    tau(ii) = 2*pi* trapz(fliplr(altitude), fliplr(re_meters.^2 .* total_Nc_meters));

end

% add a zero at the begining!
tau = [0,tau];


% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];

figure;
plot(fliplr(vr.re(index_cloudBottom:index_cloudTop2)),tau,'k'); 
set(gca,'YDir','reverse')
ylabel('$$\tau$$','interpreter','latex'); xlabel('$$r_{e}$$ $$(\mu m)$$','Interpreter','latex')
title('Cloud Top = $$max(LWC)$$', 'Interpreter','latex')
grid on; grid minor; hold on; 

% Plot the modis droplet estimate as a constant vertical line
xl0 = xline(re_17_modis,'--',['MODIS $$r_{2.1} = $$',num2str(re_17_modis), '$$\mu m$$'], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(tau_17_modis,'--',['MODIS $$\tau_{2.1} = $$',num2str(tau_17_modis)], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';


% Let's compute the mean number concentration within this cloud and print
% it on our plot
mean_Nc = mean(vr.total_Nc(index_cloudBottom:index_cloudTop));

dim = [.2 .5 .3 .3];
str = ['$$< N_c > = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18,'FontWeight','bold');



%% Compute the Liquid water path using the threshold defined by Painemal and Zudema

% Let's find all times when lwc is greater than the threshold above
% This will define each individual cloud in our time series

index_lwc = vr.lwc>lwc_lim;

index_minLength = 10;                            % each segment has to have atleast 3 data points to be counted as a cloud

% We want to find where these logical ones start and end. A 1 next to a 0
% will have a different of 1. A 1 next to a 1 or a 0 next to a 0 will have
% a difference of 0
% We add zeros to either end of our vector so that the index is the same as
% the original vector. The difference vector shrinks the length of the
% input by 1

% So the odd entries are the start of the cloud, the even entires are the
% end of the cloud

b = find(diff([0 index_lwc 0]));
c = b(2:2:end) - b(1:2:end);

d = find(c>=index_minLength);

cloud_index = cell(length(d),1);

for ii = 1:length(d)
    cloud_index{ii} = b(2*d(ii) - 1):1:b(2*d(ii))-1;
end



%% Lets plot the time series of all clouds found in the above section

figure;

colororder('default')

for ii = 1:length(cloud_index)
    hold on
    plot(vr.time(cloud_index{ii}), vr.re(cloud_index{ii}),'Color',[0, 0.4470, 0.7410]); 
   
end
xlabel('Time (sec)'); ylabel('r_{e} (\mum)')

yyaxis right

for ii = 1:length(cloud_index)
    hold on; 
    plot(vr.time(cloud_index{ii}), vr.altitude(cloud_index{ii}), 'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-'); 
end
ylabel('Altitude (m agl)')
set(gca,'ycolor',[0.8500, 0.3250, 0.0980]) 

grid on; grid minor; 

set(gcf, 'Position', [0 0 1300 400])




%% Make histogram of the size distribution at time = 6256 when the number concentration is constant


figure;
histogram('BinEdges',vr.drop_radius_bin_edges ,'BinCounts',vr.Nc(:,pz_secSinceStart)); 
xlabel('Droplet Radius (\mum)'); ylabel('N_c (cm^{-3})')
grid on; grid minor; hold on;
xlim([0,20])




%% Project the location of VOCALS-REx using a heading and a linear distance

% % Define the folder where all VOCALS-REx data is stored
% foldername = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
%     'VOCALS_REx/vocals_rex_data/SPS_1/'];
% 

% % ----- November 9 data -----
% filename = 'RF11.20081109.125700_213600.PNI.nc';
% 
% % Load data
% vocalsRex = readVocalsRex([foldername,filename]);
% 
% % find vertical profiles
% LWC_threshold = 0.03;       % g/m^3
% stop_at_max_LWC = false;
% Nc_threshold = 10;           % cm^{-3}
% vert_profs = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_LWC, Nc_threshold);

% distance travelled by the plane. We assume it's an arc length
d0 = 1000;          % m

% head realtive to due North
heading = 0;            % degrees

% set up the initial coordinates
lat0 = 0; 
long0 = 0;

% Create a World Geodetic System of 1984 (WGS84) reference ellipsoid with a length unit of meters.
wgs84 = wgs84Ellipsoid("m");

% compute the linear distance between two points on an ellipsoid
dLat = 0.00001;
Lat2Add = dLat;

d = distance(lat0, long0, lat0 + Lat2Add, long0,wgs84);

while d<1000
    
    Lat2Add = Lat2Add + dLat;
    d = distance(lat0, long0, lat0 + Lat2Add, long0,wgs84);

end


% So, a change in the latitude of 0.009 degrees results in a change of 1
% linear km

d_Lat_1km = 0.009;      % degrees to change latitude for a 1km linear change

% Assuming the Earth is a sphere...
% We can calculate the change in longitutde for a given latitude that
% results in a 1km linear change
con = physical_constants();

d_Long_1km = d0/(con.R_earth * cosd(lat0));         % degrees



% So I start with some lat and long position. To move along a sphere, I
% need to invoke spherical trigonometry. But if I incorrectly use circular
% trig, I can solve this right now. 
% 
% My justification for ignoring spherical trigonometry is that the
% circumference of the Earth is much much larger than the arc lengths
% travelled by a cloud over 20 minutes. The circumference of the Earth is
% 4e7 meters long. The clouds measured by VOCALS-REx travel an average of
% 7.5 m/s. Over 30 minutes, this would travel a distance of 13.5 km or 
% 1.3e4 meters, roughly 1/1000 th of Earth's circumference


