%% Load some VOCALS-REx Data

% Define the folder where all VOCALS-REx data is stored
foldername = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
    'VOCALS_REx/vocals_rex_data/SPS_1/'];


% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% Oct-25-2008 Data
% filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';





% Load data
vocalsRex = readVocalsRex([foldername,filename]);

% find vertical profiles
LWC_threshold = 0.005;       % g/m^3
stop_at_max_LWC = false;
Nc_threshold = 1;           % cm^{-3}
vert_profs = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_LWC, Nc_threshold);

%% Plot the VOCALS-REx flight path

% Plot MODIS measured relfectance at 650nm with VOCALS flight path
figure; geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10,'.')
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
title('VOCALS-REx flight path','Interpreter','latex', 'FontSize', 40)
geolimits([min(vocalsRex.latitude)-5 max(vocalsRex.latitude)+5],[min(vocalsRex.longitude)-5 max(vocalsRex.longitude)+5])



%% Plot total number concentration versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

figure;
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot effective radius versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure;
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re);
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$r_e$  $(\mu m)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot total number concentration versus time and altitude versus time for entire data set


figure;
semilogy((vocalsRex.time), vocalsRex.total_Nc);
grid on; grid minor;
xlabel('Time (sec since takeoff)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot((vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])


%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time
% ----- VOCALS TIME -------



figure; subplot(2,1,1)
semilogy(double(vocalsRex.time), vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(double(vocalsRex.time), vocalsRex.re);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right;
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and LWC versus time on two subplots
% Overlay the airplane altitude on both


UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format


figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.lwc);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$LWC$ $(g/m^{3})$','Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1275, 600])



%% Find data for when the plane ascends or descends cleanly through a cloud deck

% dz/dt must be non-zero.
% Total Nc has to start at a value below 1
% Total Nc has to end at a value below 1

% First lets crop the data only to those portions where the plane is
% ascending or descending

dz_dt = diff(vocalsRex.altitude)./diff(double(vocalsRex.time))';

% Compute the mean with a sliding window for every 10 data points
% This will smooth out the data and make the horizontal flight segments
% easier to find

dz_dt_mean = movmean(dz_dt,20);

% Plot the derivative
UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

time_UTC = UTC_starttime + double(vocalsRex.time)./3600;            % Time in UTC


figure;
plot(time_UTC(2:end), dz_dt)
hold on;
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])




% Link axes with total number concentration
figure; subplot(2,1,1)
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

% Plot total number concentration versus time and altitude
subplot(2,1,2)
semilogy(time_UTC, vocalsRex.total_Nc);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('Total $N_c$ $(cm^{-3})$','Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 450])


%% Plot the individual vertical profiles versus their original time

[vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex);

% Let's see what they look like
figure;
for ii=1:length(vert_profs.Nc)
    semilogy(vert_profs.time_UTC{ii}, vert_profs.Nc{ii}, '.')
    hold on;
end
grid on; grid minor
xlabel('Time (UTC)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('Vertical profiles only', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])




%% Plot the individual vertical profiles normalized by cloud top height
% Place them all on the same panel


% Let's see what they look like
figure;


% First plot the LWC
subplot(1,3,1)

for ii=1:length(vert_profs.lwc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./(max(vert_profs.altitude{ii})- min(vert_profs.altitude{ii}));

    l = semilogx(vert_profs.lwc{ii}, z_norm, 'LineWidth',3);
    % Set the transparency to 50%
    l.Color(4) = 0.5;

    hold on;

end

grid on; grid minor
ylabel('$Z/Z_{top}$','Interpreter','latex')
xlabel('$LWC$  $(g/m^{3})$','Interpreter','latex')


% Plot the effective radius next
subplot(1,3,2)

for ii=1:length(vert_profs.re)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.re{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$r_e$  $(\mu m)$','Interpreter','latex')


% Plot the total droplet number concentration last
subplot(1,3,3)

for ii=1:length(vert_profs.Nc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.Nc{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$N_c$  $(cm^{-3})$','Interpreter','latex')



set(gcf, 'Position',[0 0 1300, 450])


%% Get all indices for brushed data

% brushedData is the data grabbed from the plot

index_v = [];

for ii=1:size(brushedData,1)

    is_true = brushedData(ii,1)==vocalsRex.time;

    index_v = [index_v, is_true];

end

% sum all columns up into one logical vector!
index_v = logical(sum(index_v,2));


% Save all new data

% create structure with just vertical profiles
vert_profs.time = vocalsRex.time(index_v);
vert_profs.total_Nc = vocalsRex.total_Nc(index_v)';
vert_profs.lwc = vocalsRex.lwc(index_v)';
vert_profs.latitude = vocalsRex.latitude(index_v);
vert_profs.longitude = vocalsRex.longitude(index_v);
vert_profs.altitude = vocalsRex.altitude(index_v)';
vert_profs.re = vocalsRex.re(index_v)';


%% Plot just the vertical profiles taken from above
% Plot their LWC versus altitude


figure;
plot(vert_profs.lwc, vert_profs.altitude/1e3, '.');
grid on; grid minor;
xlabel('LWC $(g/m^3)$','Interpreter','latex')
ylabel('Altitude $(km)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])


%% Plot droplet distribution at the base of some vert profile

% Plot one distribuiton at the bottom of the cloud
index_time1 = vert_profs.time{1}(1);

% plot another distribution at the top of the cloud
index_time2 = vert_profs.time{1}(end);


figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1));
h.FaceColor = mySavedColors(11, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
    sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');
xline(re1,'k--', 'LineWidth',2)




% plot the other distribution

figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2));
h.FaceColor = mySavedColors(14, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
    sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');
xline(re2,'k--', 'LineWidth',2)





%% Plot histogram of two droplet distributions on the same figure!

% Define the min and max radius values to plot
r_min = 0;      % microns
r_max = 30;     % microns

% vert_profile to plot
index2plot = 1;

% plot a distribution at the top of the cloud
index_time1 = vert_profs.time{index2plot}(10)+1;


% plot another distribution at the top of the cloud
index_time2 = vert_profs.time{index2plot}(40)+1;


% Compute the effective radius for the two distributions and plot it as a solid vertical line
re1 = vert_profs.re{index2plot}(1);
re2 = vert_profs.re{index2plot}(end);

% re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
%       sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');
% re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
%       sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');


f1 = figure;

% Plot the distribution at cloud bottom first
h1 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1));
h1.FaceColor = mySavedColors(11, 'fixed');
h1.FaceAlpha = 0.7;
h1.EdgeAlpha = 1;
hold on
xline(re1,'k--', 'LineWidth',4)

h2 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2));
h2.FaceColor = mySavedColors(14, 'fixed');
h2.FaceAlpha = 0.4;
h2.EdgeAlpha = 1;
xline(re2,'k--', 'LineWidth',4)

xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
ylabel('$n(r)$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize',32);
grid on; grid minor; hold on;
xlim([r_min, r_max])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
set(gcf, 'Position',[0 0 1000, 600])

% Label first xline
annotation(f1,'textbox',[0.209125 0.917777777777778 0.0911875000000002 0.0711111111111111],...
    'String',['$r_e = ', num2str(re1), ' \mu m$'],...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label second x line
annotation(f1,'textbox',[0.397406250000002 0.917777777777778 0.0911875000000001 0.0711111111111111],...
    'String',['$r_e = ', num2str(re2), ' \mu m$'],...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');


legend('Cloud Bottom', '', 'Cloud Top', 'Location','best')


%% Plot droplet number distribution versus altitude as a imagesc plot

profile_idx = 1;

% imagesc assumes linearly spaced vectors along the x and y axes. Use p
% color for non-linear sampling

% plot the altitude on the y-axis
Y_data = vocalsRex.altitude(vert_profs.time{profile_idx});
% plot the droplet sizes on the x-axis
X_data = vocalsRex.drop_radius_bin_center;

[X,Y] = meshgrid(X_data, Y_data);

figure;
s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{profile_idx})');

% reduce the transparancy of the bin edges
s.EdgeAlpha = 0;

xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('Altitude ($m$)', 'Interpreter','latex')
c = colorbar;

xlim([0 50])
set(gca, 'YDir', 'normal')

ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');

set(gcf, 'Position', [0 0 800 600])

set(gca,'ColorScale','log')

set(gca,'CLim', [10^(-1) , 10^2])


%% Plot droplet number distribution versus altitude as a imagesc plot - SMOOTHED VERSION

profile_idx = 1;


interp_r = linspace(vocalsRex.drop_radius_bin_center(1), vocalsRex.drop_radius_bin_center(end), 10000);

d = zeros(size(vocalsRex.Nc(:, vert_profs.time{profile_idx})',1), numel(interp_r));

for ii = 1:size(d,1)
    d(ii,:) = interp1(vocalsRex.drop_radius_bin_center, vocalsRex.Nc(:,vert_profs.time{profile_idx}(ii))', interp_r);
end


% plot the altitude on the y-axis
Y_data = vocalsRex.altitude(vert_profs.time{profile_idx});
% plot the droplet sizes on the x-axis
X_data = interp_r;

[X,Y] = meshgrid(X_data, Y_data);

figure;
s = pcolor(X,Y,d);

% reduce the transparancy of the bin edges
s.EdgeAlpha = 0;


xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('Altitude ($m$)', 'Interpreter','latex')
c = colorbar;

xlim([0 50])
set(gca, 'YDir', 'normal')

ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');

set(gcf, 'Position', [0 0 800 600])

% Set the colorscale to be logarithmic
set(gca,'ColorScale','log')

set(gca,'CLim', [10^(-1) , 10^2])



%% Fit a distribution to the droplet distribution data

% set a minimum value so the there are no 0's in my data
min_nr_value = 10^(-3);

r = double(vocalsRex.drop_radius_bin_center);

index2plot = 1;
time2plot = 10;

% read in the number concentration data
% add 1 index when reading the Vocals rex data because the vocalsRex time
% starts at 0
Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot}(time2plot)+1);     % cm^(-3) - total number of droplets between some radius range

nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

r = r(nr_data>min_nr_value);
nr_data = nr_data(nr_data>min_nr_value);

% lets find the modal radius
[~,idx] = max(nr_data);
r_modal = r(idx);

% define the width of the distribution.
% this value has to be an integer >= 1
std_dev = 26;

% Compute the total the number concentration
Nc = vert_profs.Nc{index2plot}(time2plot);

% compute the size distribution
[nr1, r1] = gamma_size_distribution_kokhanovsky(r_modal, std_dev, Nc);

% perform a gamma distribution fit
[nr_fit, mu_fit, rms_min] = fit_gamma_size_distribution_kokhanovsky(nr_data, r);


figure; plot(r, nr_data)
hold on
plot(r1, nr1)
plot(r, nr_fit)
grid on; grid minor
xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('$n(r)$ $(cm^{-3} \mu m^{-3})$', 'Interpreter','latex')
legend('VOCALS',['Manual, \mu = ', num2str(std_dev)], ['rms fit, \mu = ', num2str(mu_fit)])


%% Testing gamma size distribution fits using Matlab

min_nr_value = 10^(-3);

r = double(vocalsRex.drop_radius_bin_center)';

index2plot = 1;
time2plot = 10;

% read in the number concentration data
% add 1 index when reading the Vocals rex data because the vocalsRex time
% starts at 0
Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot}(time2plot)+1);     % cm^(-3) - total number of droplets between some radius range

nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

r = r(nr_data>min_nr_value);
nr_data = nr_data(nr_data>min_nr_value);

% lets find the modal radius
[~,idx] = max(nr_data);
r_modal = r(idx);


% Compute the total the number concentration
Nc = vert_profs.Nc{index2plot}(time2plot);       % cm^(-3)

%gamma_fit = fitdist(r, 'gamma', 'Frequency', nr_data);
% logNormal_fit = fitdist(r, 'Lognormal', 'Frequency', nr_data);

gamma_fit = fitdist(nr_data, 'gamma');
% logNormal_fit = fitdist(nr_data, 'Lognormal');

%gamma_fit = makedist("Gamma","a",r_modal, "b", sigma0);

% plot the data and all the distributions

f = figure;
p1 = plot(r, nr_data);
xlabel('r')
ylabel('n(r)')
grid on; grid minor

hold on
plot(r, Nc*pdf(gamma_fit, r));

%plot(r, pdf(logNormal_fit, r));

legend('data','gamma fit','log normal fit')





%% Testing Matlabs built-in gamma size distributions

r_modal = 5;           % microns

% define the effective variance
sigma0 = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(sigma0)

    gamma_fit = makedist("Gamma","a",r_modal, "b", sigma0(nn));
    plot(r, Nc*pdf(gamma_fit, r));
    hold on
    legend_str{nn} = ['\sigma_0 = ', num2str(sigma0(nn))];
end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Matlabs Gamma distribution')
legend(legend_str)


%% testing my own gamma size distribution using libRadTrans equation


r_modal = 5;           % microns

% define the effective variance
sigma0 = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(sigma0)

    [nr1, r1] = gamma_size_distribution_libRadTran(5, sigma0(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\sigma_0 = ', num2str(sigma0(nn))];

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('LibRadTran Gamma distribution')
legend(legend_str)

%% testing my own gamma size distribution using Kokhanovsky's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 1:10;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = gamma_size_distribution_kokhanovsky(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Kokhanovsky Gamma distribution')
legend(legend_str)


%% testing my own gamma size distribution using Pruppacher's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)

% test the integral - should equal Nc
total_Nc = zeros(1, numel(std_dev));

% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = gamma_size_distribution_pruppacher(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];
    total_Nc(nn) = trapz(r1, nr1);

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Pruppacher Gamma distribution')
legend(legend_str)



%% testing my own log noraml size distribution using Kokhanovsky's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)

% test the integral - should equal Nc
total_Nc = zeros(1, numel(std_dev));

% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = lognormal_size_distribution_kokhanovsky(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];
    total_Nc(nn) = trapz(r1, nr1);

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Kokhanovsky Lognormal distribution')
legend(legend_str)




%% Load some VOCALS-REx Data

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% November 11, 2008 data
%filename = 'RF12.20081111.125000_214500.PNI.nc';


% ----- November 9 data -----
foldername = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/2008_11_09/';
filename = 'RF11.20081109.125700_213600.PNI.nc';



% Load data
vocalsRex = readVocalsRex([foldername,filename]);

%% Plot the VOCALS-REx flight path

% Plot MODIS measured relfectance at 650nm with VOCALS flight path
figure; geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10,'.')
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
title('VOCALS-REx flight path','Interpreter','latex', 'FontSize', 40)
geolimits([min(vocalsRex.latitude)-5 max(vocalsRex.latitude)+5],[min(vocalsRex.longitude)-5 max(vocalsRex.longitude)+5])



%% Plot total number concentration versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

figure;
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot effective radius versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure;
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re);
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$r_e$  $(\mu m)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot total number concentration versus time and altitude versus time for entire data set


figure;
semilogy((vocalsRex.time), vocalsRex.total_Nc);
grid on; grid minor;
xlabel('Time (sec since takeoff)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot((vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])


%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time
% ----- VOCALS TIME -------



figure; subplot(2,1,1)
semilogy(double(vocalsRex.time), vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(double(vocalsRex.time), vocalsRex.re);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right;
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and LWC versus time on two subplots
% Overlay the airplane altitude on both


UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format


figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc);
grid on; grid minor;
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.lwc);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('$LWC$ $(g/m^{3})$','Interpreter','latex')

hold on;
yyaxis right;
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1275, 600])



%% Find data for when the plane ascends or descends cleanly through a cloud deck

% dz/dt must be non-zero.
% Total Nc has to start at a value below 1
% Total Nc has to end at a value below 1

% First lets crop the data only to those portions where the plane is
% ascending or descending

dz_dt = diff(vocalsRex.altitude)./diff(double(vocalsRex.time))';

% Compute the mean with a sliding window for every 10 data points
% This will smooth out the data and make the horizontal flight segments
% easier to find

dz_dt_mean = movmean(dz_dt,20);

% Plot the derivative
UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

time_UTC = UTC_starttime + double(vocalsRex.time)./3600;            % Time in UTC


figure;
plot(time_UTC(2:end), dz_dt)
hold on;
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])




% Link axes with total number concentration
figure; subplot(2,1,1)
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor;
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

% Plot total number concentration versus time and altitude
subplot(2,1,2)
semilogy(time_UTC, vocalsRex.total_Nc);
grid on; grid minor;
xlabel('Time (hours)','Interpreter','latex')
ylabel('Total $N_c$ $(cm^{-3})$','Interpreter','latex')

hold on;
yyaxis right;
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 450])


%% Plot the individual vertical profiles versus their original time

[vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex);

% Let's see what they look like
figure;
for ii=1:length(vert_profs.Nc)
    semilogy(vert_profs.time_UTC{ii}, vert_profs.Nc{ii}, '.')
    hold on;
end
grid on; grid minor
xlabel('Time (UTC)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('Vertical profiles only', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])




%% Plot the individual vertical profiles normalized by cloud top height
% Place them all on the same panel


% Let's see what they look like
figure;


% First plot the LWC
subplot(1,3,1)

for ii=1:length(vert_profs.lwc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./(max(vert_profs.altitude{ii})- min(vert_profs.altitude{ii}));

    l = semilogx(vert_profs.lwc{ii}, z_norm, 'LineWidth',3);
    % Set the transparency to 50%
    l.Color(4) = 0.5;

    hold on;

end

grid on; grid minor
ylabel('$Z/Z_{top}$','Interpreter','latex')
xlabel('$LWC$  $(g/m^{3})$','Interpreter','latex')


% Plot the effective radius next
subplot(1,3,2)

for ii=1:length(vert_profs.re)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.re{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$r_e$  $(\mu m)$','Interpreter','latex')


% Plot the total droplet number concentration last
subplot(1,3,3)

for ii=1:length(vert_profs.Nc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.Nc{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$N_c$  $(cm^{-3})$','Interpreter','latex')



set(gcf, 'Position',[0 0 1300, 450])


%% Get all indices for brushed data

% brushedData is the data grabbed from the plot

index_v = [];

for ii=1:size(brushedData,1)

    is_true = brushedData(ii,1)==vocalsRex.time;

    index_v = [index_v, is_true];

end

% sum all columns up into one logical vector!
index_v = logical(sum(index_v,2));


% Save all new data

% create structure with just vertical profiles
vert_profs.time = vocalsRex.time(index_v);
vert_profs.total_Nc = vocalsRex.total_Nc(index_v)';
vert_profs.lwc = vocalsRex.lwc(index_v)';
vert_profs.latitude = vocalsRex.latitude(index_v);
vert_profs.longitude = vocalsRex.longitude(index_v);
vert_profs.altitude = vocalsRex.altitude(index_v)';
vert_profs.re = vocalsRex.re(index_v)';


%% Plot just the vertical profiles taken from above
% Plot their LWC versus altitude


figure;
plot(vert_profs.lwc, vert_profs.altitude/1e3, '.');
grid on; grid minor;
xlabel('LWC $(g/m^3)$','Interpreter','latex')
ylabel('Altitude $(km)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])


%% Plot droplet distribution at the base of some vert profile

% Plot one distribuiton at the bottom of the cloud
index_time1 = vert_profs.time{1}(1);

% plot another distribution at the top of the cloud
index_time2 = vert_profs.time{1}(end);


figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1));
h.FaceColor = mySavedColors(11, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
    sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');
xline(re1,'k--', 'LineWidth',2)




% plot the other distribution

figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2));
h.FaceColor = mySavedColors(14, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
    sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');
xline(re2,'k--', 'LineWidth',2)





%% Plot histogram of two droplet distributions on the same figure!

% Define the min and max radius values to plot
r_min = 0;      % microns
r_max = 30;     % microns

% vert_profile to plot
index2plot = 1;

% plot a distribution at the top of the cloud
index_time1 = vert_profs.time{index2plot}(10)+1;


% plot another distribution at the top of the cloud
index_time2 = vert_profs.time{index2plot}(40)+1;


% Compute the effective radius for the two distributions and plot it as a solid vertical line
re1 = vert_profs.re{index2plot}(1);
re2 = vert_profs.re{index2plot}(end);

% re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
%       sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');
% re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
%       sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');


f1 = figure;

% Plot the distribution at cloud bottom first
h1 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1));
h1.FaceColor = mySavedColors(11, 'fixed');
h1.FaceAlpha = 0.7;
h1.EdgeAlpha = 1;
hold on
xline(re1,'k--', 'LineWidth',4)

h2 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2));
h2.FaceColor = mySavedColors(14, 'fixed');
h2.FaceAlpha = 0.4;
h2.EdgeAlpha = 1;
xline(re2,'k--', 'LineWidth',4)

xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
ylabel('$n(r)$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize',32);
grid on; grid minor; hold on;
xlim([r_min, r_max])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
set(gcf, 'Position',[0 0 1000, 600])

% Label first xline
annotation(f1,'textbox',[0.209125 0.917777777777778 0.0911875000000002 0.0711111111111111],...
    'String',['$r_e = ', num2str(re1), ' \mu m$'],...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label second x line
annotation(f1,'textbox',[0.397406250000002 0.917777777777778 0.0911875000000001 0.0711111111111111],...
    'String',['$r_e = ', num2str(re2), ' \mu m$'],...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');


legend('Cloud Bottom', '', 'Cloud Top', 'Location','best')


%% Plot droplet number distribution versus altitude as a imagesc plot

profile_idx = 1;

% imagesc assumes linearly spaced vectors along the x and y axes. Use p
% color for non-linear sampling

% plot the altitude on the y-axis
Y_data = vocalsRex.altitude(vert_profs.time{profile_idx});
% plot the droplet sizes on the x-axis
X_data = vocalsRex.drop_radius_bin_center;

[X,Y] = meshgrid(X_data, Y_data);

figure;
s = pcolor(X,Y,vocalsRex.Nc(:, vert_profs.time{profile_idx})');

% reduce the transparancy of the bin edges
s.EdgeAlpha = 0;

xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('Altitude ($m$)', 'Interpreter','latex')
c = colorbar;

xlim([0 50])
set(gca, 'YDir', 'normal')

ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');

set(gcf, 'Position', [0 0 800 600])

set(gca,'ColorScale','log')

set(gca,'CLim', [10^(-1) , 10^2])


%% Plot droplet number distribution versus altitude as a imagesc plot - SMOOTHED VERSION

profile_idx = 1;


interp_r = linspace(vocalsRex.drop_radius_bin_center(1), vocalsRex.drop_radius_bin_center(end), 10000);

d = zeros(size(vocalsRex.Nc(:, vert_profs.time{profile_idx})',1), numel(interp_r));

for ii = 1:size(d,1)
    d(ii,:) = interp1(vocalsRex.drop_radius_bin_center, vocalsRex.Nc(:,vert_profs.time{profile_idx}(ii))', interp_r);
end


% plot the altitude on the y-axis
Y_data = vocalsRex.altitude(vert_profs.time{profile_idx});
% plot the droplet sizes on the x-axis
X_data = interp_r;

[X,Y] = meshgrid(X_data, Y_data);

figure;
s = pcolor(X,Y,d);

% reduce the transparancy of the bin edges
s.EdgeAlpha = 0;


xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('Altitude ($m$)', 'Interpreter','latex')
c = colorbar;

xlim([0 50])
set(gca, 'YDir', 'normal')

ylabel(c,'$n(r)$ $(m^{-3}$)','FontSize',25, 'interpreter', 'latex');

set(gcf, 'Position', [0 0 800 600])

% Set the colorscale to be logarithmic
set(gca,'ColorScale','log')

set(gca,'CLim', [10^(-1) , 10^2])



%% Fit a distribution to the droplet distribution data

% set a minimum value so the there are no 0's in my data
min_nr_value = 10^(-3);

r = double(vocalsRex.drop_radius_bin_center);

index2plot = 1;
time2plot = 10;

% read in the number concentration data
% add 1 index when reading the Vocals rex data because the vocalsRex time
% starts at 0
Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot}(time2plot)+1);     % cm^(-3) - total number of droplets between some radius range

nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

r = r(nr_data>min_nr_value);
nr_data = nr_data(nr_data>min_nr_value);

% lets find the modal radius
[~,idx] = max(nr_data);
r_modal = r(idx);

% define the width of the distribution.
% this value has to be an integer >= 1
std_dev = 26;

% Compute the total the number concentration
Nc = vert_profs.Nc{index2plot}(time2plot);

% compute the size distribution
[nr1, r1] = gamma_size_distribution_kokhanovsky(r_modal, std_dev, Nc);

% perform a gamma distribution fit
[nr_gamma_fit, mu_fit, var_coef, rms_min_gamma] = fit_gamma_size_distribution_kokhanovsky(nr_data, r);

% perform a lognormal distribution fit
[nr_lognormal_fit, sigma_fit, rms_min_lognormal] = fit_lognormal_size_distribution_kokhanovsky(nr_data, r);

% show the rms values
disp(['Gamma rms fit = ',num2str(rms_min_gamma), '    Lognormal rms fit = ',num2str(rms_min_lognormal)])


figure; plot(r, nr_data)
hold on
plot(r1, nr1)
plot(r, nr_gamma_fit)
plot(r, nr_lognormal_fit)
grid on; grid minor
xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('$n(r)$ $(cm^{-3} \mu m^{-3})$', 'Interpreter','latex')
legend('VOCALS',['Manual, \mu = ', num2str(std_dev)], ['gamma fit, \mu = ', num2str(mu_fit)],...
    ['lognormal fit, \sigma = ', num2str(sigma_fit)])



%% How does the variance of the distribution change with altitude within cloud?

% set a minimum value so the there are no 0's in my data
min_nr_value = 10^(-3);

r = double(vocalsRex.drop_radius_bin_center);

index2plot = 1;
time2plot = 10;

n_time = numel(vert_profs.time{index2plot});

% read in the number concentration data
% add 1 index when reading the Vocals rex data because the vocalsRex time
% starts at 0

for nn = 1:n_time
    Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot}(nn)+1);     % cm^(-3) - total number of droplets between some radius range

    nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

    r = r(nr_data>min_nr_value);
    nr_data = nr_data(nr_data>min_nr_value);



    % perform a gamma distribution fit
    [nr_gamma_fit(nn), mu_fit(nn), var_coef(nn), rms_min_gamma(nn)] = fit_gamma_size_distribution_kokhanovsky(nr_data, r);

    % perform a lognormal distribution fit
    [nr_lognormal_fit(nn), sigma_fit(nn), rms_min_lognormal(nn)] = fit_lognormal_size_distribution_kokhanovsky(nr_data, r);

end


figure; plot(var_coef, vert_profs.altitude{index2plot})
hold on
plot(sigma_fit, vert_profs.altitude{index2plot})

grid on; grid minor
xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
ylabel('$n(r)$ $(cm^{-3} \mu m^{-3})$', 'Interpreter','latex')
legend('VOCALS',['Manual, \mu = ', num2str(std_dev)], ['gamma fit, \mu = ', num2str(mu_fit)],...
    ['lognormal fit, \sigma = ', num2str(sigma_fit)])

%% Testing gamma size distribution fits using Matlab

min_nr_value = 10^(-3);

r = double(vocalsRex.drop_radius_bin_center)';

index2plot = 1;
time2plot = 49;

% read in the number concentration data
% add 1 index when reading the Vocals rex data because the vocalsRex time
% starts at 0
Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot}(time2plot)+1);     % cm^(-3) - total number of droplets between some radius range

nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

r = r(nr_data>min_nr_value);
nr_data = nr_data(nr_data>min_nr_value);

% lets find the modal radius
[~,idx] = max(nr_data);
r_modal = r(idx);


% Compute the total the number concentration
Nc = vert_profs.Nc{index2plot}(time2plot);       % cm^(-3)

%gamma_fit = fitdist(r, 'gamma', 'Frequency', nr_data);
% logNormal_fit = fitdist(r, 'Lognormal', 'Frequency', nr_data);

%gamma_fit = fitdist(nr_data, 'gamma');
% logNormal_fit = fitdist(nr_data, 'Lognormal');

%gamma_fit = makedist("Gamma","a",r_modal, "b", sigma0);

% make a gamma distribution using Kokhavonvsky
[nr_gamma, r_gamma] = gamma_size_distribution_kokhanovsky(r_modal, 110, Nc);

[nr_fit, fit_params] = fit_gamma_size_distribution_kokhanovsky(nr_data, r);

% plot the data and all the distributions

f = figure;
p1 = plot(r, nr_data);
xlabel('r')
ylabel('n(r)')
grid on; grid minor

hold on
xline(vert_profs.re{index2plot}(time2plot), 'k--', 'Label','r_e','LineWidth',3, ...
    'FontSize',15)

%plot(r, Nc*pdf(gamma_fit, r));

%plot(r, pdf(logNormal_fit, r));

% plot the manual gamma distribution fit
plot(r_gamma, nr_gamma);

% plot the algorithmic gamma distribution fit
plot(r, nr_fit);

legend('data', 'r_e', 'gamma fit - manual', 'gamma fit - algorithm')





%% Testing Matlabs built-in gamma size distributions

r_modal = 5;           % microns

% define the effective variance
sigma0 = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(sigma0)

    gamma_fit = makedist("Gamma","a",r_modal, "b", sigma0(nn));
    plot(r, Nc*pdf(gamma_fit, r));
    hold on
    legend_str{nn} = ['\sigma_0 = ', num2str(sigma0(nn))];
end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Matlabs Gamma distribution')
legend(legend_str)


%% testing my own gamma size distribution using libRadTrans equation


r_modal = 5;           % microns

% define the effective variance
sigma0 = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(sigma0)

    [nr1, r1] = gamma_size_distribution_libRadTran(5, sigma0(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\sigma_0 = ', num2str(sigma0(nn))];

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('LibRadTran Gamma distribution')
legend(legend_str)

%% testing my own gamma size distribution using Kokhanovsky's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 1:10;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = gamma_size_distribution_kokhanovsky(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Kokhanovsky Gamma distribution')
legend(legend_str)


%% testing my own gamma size distribution using Pruppacher's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)

% test the integral - should equal Nc
total_Nc = zeros(1, numel(std_dev));

% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = gamma_size_distribution_pruppacher(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];
    total_Nc(nn) = trapz(r1, nr1);

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Pruppacher Gamma distribution')
legend(legend_str)



%% testing my own log noraml size distribution using NIST's equation


r_modal = 5;           % microns

% define the effective variance
std_dev = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)

% test the integral - should equal Nc
total_Nc = zeros(1, numel(std_dev));

% plot the data and all the distributions

f = figure;
for nn = 1:numel(std_dev)

    [nr1, r1] = lognormal_size_distribution_kokhanovsky(5, std_dev(nn), 100);
    plot(r1, nr1)
    hold on
    legend_str{nn} = ['\mu = ', num2str(std_dev(nn))];
    total_Nc(nn) = trapz(r1, nr1);

end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Kokhanovsky Lognormal distribution')
legend(legend_str)


%% Testing Matlabs built-in log normal size distributions

r_modal = 5;           % microns

% define the effective variance
sigma0 = 0.1:0.1:1;

% Compute the total the number concentration
Nc = 100;       % cm^(-3)


% plot the data and all the distributions

f = figure;
for nn = 1:numel(sigma0)

    lognormal_fit = makedist("Lognormal","mu",log(r_modal), "sigma", sigma0(nn));
    plot(r, Nc*pdf(lognormal_fit, r));
    hold on
    legend_str{nn} = ['\sigma_0 = ', num2str(sigma0(nn))];
    disp(['meadian = ', num2str(lognormal_fit.median)])
end

xlabel('r')
ylabel('n(r)')
grid on; grid minor
title('Matlabs Lognormal distribution')
legend(legend_str)


%% Plot the modal radius versus the effective radius


min_nr_value = 10^(-3);

r0 = double(vocalsRex.drop_radius_bin_center)';

index2plot = 3;

clear fit_parameters

for nn = 1:numel(index2plot)

    for mm = 1:numel(vert_profs.time{index2plot(nn)})

        % read in the number concentration data
        % add 1 index when reading the Vocals rex data because the vocalsRex time
        % starts at 0
        Nc_data = vocalsRex.Nc(:, vert_profs.time{index2plot(nn)}(mm)+1);     % cm^(-3) - total number of droplets between some radius range

        nr_data = Nc_data./diff(vocalsRex.drop_radius_bin_edges)';               % cm^(-3) * microns^(-1) - droplet distribution data

        r = r0(nr_data>min_nr_value);
        nr_data = nr_data(nr_data>min_nr_value);

        [~, fit_parameters(nn,mm)] = fit_gamma_size_distribution_kokhanovsky(nr_data, r);

    end
end


% plot the results
figure; subplot(1,2,1)
for ii=1:numel(vert_profs.time{index2plot})
    plot(fit_parameters(ii).r_modal./vert_profs.re{index2plot}(ii), vert_profs.altitude{index2plot}(ii) ,'.')
    hold on
end
xlabel('r_{modal}/r_e')
ylabel('Altitude (m)')
grid on; grid minor
title(['Vertical Profile ', num2str(index2plot)], 'Interpreter', 'latex')


subplot(1,2,2)
for ii=1:numel(vert_profs.time{index2plot})
    plot(fit_parameters(ii).mu, vert_profs.altitude{index2plot}(ii) ,'.')
    hold on
end
xlabel('\mu (variance of sorts)')
ylabel('Altitude (m)')
grid on; grid minor

set(gcf, 'Position',[0 0 1000 500])


%% plot a vertical profile histogram and a distribution fit

% Define the min and max radius values to plot
r_min = 0;      % microns
r_max = 30;     % microns

% define minimum number of droplet density
min_density = 10^(-1);
max_density = 10^2;

index2plot = 1;                                         % plot the first profile
time2plot = size(vert_profs.nr{index2plot}, 2);         % plot the cloud top 

% define the data used for fitting
data2plot = vert_profs.nr{index2plot}(:,time2plot);
% get rid of values below the minimum density threshold
idx2plot = data2plot>min_density;
data2plot = data2plot(idx2plot);

independent_data2plot = vert_profs.drop_radius_bin_center;
% only keep the same values from above
independent_data2plot = independent_data2plot(idx2plot);
% grab the bin edges to plot
idx2plot_num = find(idx2plot);
binEdges2plot = vert_profs.drop_radius_bin_edges([idx2plot_num; idx2plot_num(end)+1]);      

% compute the total number of droplets
N0 = vert_profs.Nc{index2plot}(time2plot);





% Fit a log noraml distribution using my function
[log_norm_fit, fit_params] = fit_lognormal_size_distribution_kokhanovsky(data2plot,independent_data2plot);


% fit a log normal distribution using Matlabs built in function
logNormal_fit = fitdist(data2plot, 'Lognormal');

% fit a normal distribution using Matlabs built in function
% first, lets define the gaussian equation
normEqn = 'a*exp(-((x-b)/c)^2)+d';
% define the first guesses for each variable in the above equation
[max_val, idx_maxVal] = max(data2plot);
start_points = [max_val, independent_data2plot(idx_maxVal), 1, 0];
[normal_fit, gof] = fit(independent_data2plot', data2plot, normEqn, 'StartPoint', start_points);

% plot the results

figure;
% Plot the distribution at cloud bottom first
h1 = histogram('BinEdges',binEdges2plot ,'BinCounts',data2plot);
h1.FaceColor = mySavedColors(1, 'fixed');
h1.FaceAlpha = 0.7;
h1.EdgeAlpha = 1;
hold on

% what profile are we plotting?
title(['Vertical Profile ', num2str(index2plot)], 'Interpreter','latex')

% set axes limits and labels
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
ylabel('$n(r)$ ($cm^{-3} \, \mu m^{-1}$)', 'Interpreter','latex', 'FontSize',32);
grid on; grid minor; hold on;
xlim([r_min, r_max])
ylim([min_density max_density])
%set(gca, 'YScale', 'log')
set(gcf, 'Position',[0 0 1000, 600])


% Plot the log noraml distribution fit using my function
hold on
plot(independent_data2plot, log_norm_fit)

% plot Matlab's lognormal distribution fit  
plot(independent_data2plot, pdf(logNormal_fit, independent_data2plot));


% plot a random lognormal distribution 
hold on
[log_norm_dist, r] = lognormal_size_distribution_kokhanovsky(8, 0.2, N0);
plot(r, log_norm_dist)


% plot a normal distribution 
hold on
plot(independent_data2plot(1):0.01:independent_data2plot(end), normal_fit(independent_data2plot(1):0.01:independent_data2plot(end)))

legend('n(r)', 'my lognoraml fit','MATLAB lognoraml fit','lognormal dist', 'normal dist')
