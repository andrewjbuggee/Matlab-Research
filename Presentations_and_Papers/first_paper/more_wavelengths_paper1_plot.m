% Contour plot for paper 1

% I can't think of a better way to look at this. Peter wants me to include
% all three variables. His view of using all three variables by taking the
% difference between r_top and r_bot didn't have the same effect as when I
% a 2D slice of the rms at the minimum tau value. 

% What is it I'm trying to show? That when the number of wavelengths used
% in the retrieval are increased, the space of possible solutions shrinks,
% especially along the r_bot dimension. How do I define the space of
% possible solutions? It should be the region where the rms difference
% between the measurements and the forward model calculations are less than
% the rms of the measurement uncertainty. 

% But let's rethink the solution space. How is it defined? Yes, I could set
% the threshold defined above in my iterative solver, but I haven't. How I
% currently defined convergence? 

% To Peter's point, you can't expect to do that much better than the
% uncertainty of the measurements. It doesn't make sense to do so because
% within the uncertainty range, you don't know what is true and what isn't.

% Do I need the number values on my contour plot? Isn't it more important
% to show there is a bullseye? This only works if there isn't a bullseye
% when using 7 MODIS wavelengths. But there will be. What is the slope of
% the bullseye is steeper when using more wavelengths? Look at the rms
% surface as a 3D plot where z is the rms value, x and y are the radii are
% r-top and r-bot


% By Andrew John Buggee

clear variables


%% LOAD DATA SET

%load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-26-Nov-2024_rev1.mat')
%load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-27-Nov-2024_rev1.mat')
load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-11-Dec-2024_rev2.mat')

%% Want to use real EMIT geometry inputs?

% Load modis data and create input structure

% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';



    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    emitPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder2save = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/reflectance_uniqueness/'];

elseif strcmp(whatComputer,'curc')==true

    % ------ Folders on the CU Supercomputer /projects folder --------

    % Define the EMIT data folder path

    emitPath = '/projects/anbu8374/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder2save = '/scratch/alpine/anbu8374/EMIT_reflectance_uniqueness/';

    if ~exist(folder2save, 'dir')

            mkdir(folder2save)
    end





end

%%
% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitFolder = '17_Jan_2024_coast/';
%emitFolder = '17_Jan_2024_ocean/';


[emit,L1B_fileName] = retrieveEMIT_data([emitPath, emitFolder]);


% Define an index to use
%modis_idx = 110292;     % for 9 nov 2008

% 17_Jan_2024_coast - large optical depth
% row = 1112;
% col = 974;

% 17_Jan_2024_coast - small optical depth
% row = 912;
% col = 929;

% 17_Jan_2024_coast - my TBLUT algorithm found an optical depth of 6.6
pixels2use.row = 932;
pixels2use.col = 960;

% 17_Jan_2024_coast - optical depth of 3.2 and 3.8
% pixels2use.row = [932];
% pixels2use.col = [970];

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));


%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, pixels2use);

%% compute the rms of the EMIT reflectance uncertainty

Refl_emit_uncertainty = 0.05 .* Refl_emit;
rms_uncert = sqrt(mean(Refl_emit_uncertainty.^2));

% compute the rms of the EMIT reflectance uncertainty using frist 7 MODIS
% wavelenghts
wl_MODIS7_idx = [1, 4, 6, 7, 19, 23, 29];
rms_uncert_MODIS7 = sqrt(mean(Refl_emit_uncertainty(wl_MODIS7_idx).^2));

%% Find the states with the lowest rms residul

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);


% find n smallest rms states
n_states = 50;

% store the state values at each minimum
r_top_min = zeros(n_states, 1);
r_bot_min = zeros(n_states, 1);
tau_c_min = zeros(n_states, 1);

% store the rms value and the index
min_val = zeros(n_states, 1);
idx_min = zeros(n_states, 1);


% create a new array where the rms_residual can be used to determine the
% smallest values. We have to insert a nan each time
rms_residual_placeHolder = rms_residual;


for nn = 1:n_states

    % find the smallest rms residual value, omitting nans
    [min_val(nn), idx_min(nn)] = min(rms_residual_placeHolder, [], 'all', 'omitnan');

    r_top_min(nn) = R_top_fine(idx_min(nn));
    r_bot_min(nn) = R_bot_fine(idx_min(nn));
    tau_c_min(nn) = Tau_c_fine(idx_min(nn));

    % set the minimum value to nan and omit
    rms_residual_placeHolder(idx_min(nn)) = nan;


end


% save the reflectance estimates associated with the minimum rms value
% across all three variables
min_Refl_model_fine = reshape(Refl_model_fine(r_top_fine==r_top_min(1), r_bot_fine==r_bot_min(1), tau_c_fine==tau_c_min(1),:), [],1);


%% Find the states with the lowest rms residul using the first 7 MODIS wavelengths

% Let's now seperate out the interpolated relfectance at the seven MODIS
% wavelengths

Refl_model_fine_MODIS7 = Refl_model_fine(:,:,:, wl_MODIS7_idx);


% find n smallest rms states
n_states = 50;

% store the state values at each minimum
r_top_min_MODIS7 = zeros(n_states, 1);
r_bot_min_MODIS7 = zeros(n_states, 1);
tau_c_min_MODIS7 = zeros(n_states, 1);

% store the rms value and the index
min_val_MODIS7 = zeros(n_states, 1);
idx_min_MODIS7 = zeros(n_states, 1);


% create a new array where the rms_residual can be used to determine the
% smallest values. We have to insert a nan each time
rms_residual_placeHolder_MODIS7 = rms_residual_MODIS7;


for nn = 1:n_states

    % find the smallest rms residual value, omitting nans
    [min_val_MODIS7(nn), idx_min_MODIS7(nn)] = min(rms_residual_placeHolder_MODIS7, [], 'all', 'omitnan');

    r_top_min_MODIS7(nn) = R_top_fine(idx_min_MODIS7(nn));
    r_bot_min_MODIS7(nn) = R_bot_fine(idx_min_MODIS7(nn));
    tau_c_min_MODIS7(nn) = Tau_c_fine(idx_min_MODIS7(nn));

    % set the minimum value to nan and omit
    rms_residual_placeHolder_MODIS7(idx_min_MODIS7(nn)) = nan;


end


% save the reflectance estimates associated with the minimum rms value
% across all three variables
min_Refl_model_fine_MODIS7 = reshape(Refl_model_fine_MODIS7(r_top_fine==r_top_min_MODIS7(1),...
    r_bot_fine==r_bot_min_MODIS7(1), tau_c_fine==tau_c_min_MODIS7(1),:), [],1);



%%





% -------------------------------------
% ----------- EMIT PLOTS -------------
% -------------------------------------








%% Plot 50 points with the lowest rms value

figure; plot3(r_bot_min, r_top_min, tau_c_min, '.');
% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title('50 points with lowest rms', 'FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

grid on; grid minor

% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
set(gcf, 'Position', [0 0 900 900])




%% Make 3D plot of slices along the tau dimension




% plot the boundaries of the data cube
% r_bot_slice = [max(r_bot_fine), r_bot_min(1)];
r_bot_slice = [max(r_bot_fine), ceil(r_bot_fine(end)/2)];
r_top_slice = [max(r_top_fine)];

% define the slices along the tau dimension to plot
tau_slice = [min(tau_c_fine), tau_c_min(1)];

f = figure; 

s = slice(R_bot_fine, R_top_fine, Tau_c_fine, rms_residual./rms_uncert, r_bot_slice, r_top_slice, tau_slice);

cb = colorbar('Position',[0.89108 0.11 0.016 0.815]);
% create colorbar label
ylabel(cb, '$RMS(R(\vec{x}) - \vec{m})/RMS(\delta \vec{m})$', 'FontSize', 30, 'Interpreter', 'latex')
clim([0, 7])

% create colorbar label
%ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')

% set the edge alpha to a value near 0
for nn = 1:length(s)
    s(nn).EdgeAlpha = 0.2;
end



% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title('RMS Residual using 35 EMIT Spectral Measurements','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Set position and size of figure
set(gcf, 'Position', [0 0 900 900])


%% Make 3D plot of contour slices along the tau dimension



% rms residual values to plot
lvls = [0, 0.25, 0.5:0.25:4];

% define the slices along the tau dimension to plot
tau_slice = [6, 6.2, tau_c_min(1), 6.6];

% Create figure
figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2,'Position',[0.11068 0.11 0.7254 0.815]);
hold(axes1,'on');


s = contourslice(R_bot_fine, R_top_fine, Tau_c_fine, rms_residual./rms_uncert, [], [], tau_slice, lvls);
view(3)
grid on; grid minor

hold on; plot3(r_bot_min(1), r_top_min(1), tau_c_min(1), 'r.', 'markersize', 20)

% Create colorbar
cb = colorbar(axes1,'Position',[0.89108 0.11 0.016 0.815]);
% create colorbar label
ylabel(cb, '$RMS(R(\vec{x}) - \vec{m})/RMS(\delta \vec{m})$', 'FontSize', 30, 'Interpreter', 'latex')

% set all contour lines to a certain thickness
for nn = 1:length(s)
    s(nn).LineWidth = 4;
end


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title('Using 35 EMIT Spectral Measurements','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Set position and size of figure
set(gcf, 'Position', [0 0 1000 1000])

% set 3D point of view and grid lines
view(axes1,[-56.6186981596273 23.687116207586]);
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on');



%% Create filled Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
% plot the RMS residual at the minimum optical depth and let the radii at
% cloud top and bottom varry


% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC = tau_c_fine == tau_c_min(1);
%idx_tauC = tau_c_fine == 6.2;

% Create figure
figure;


% Create axes
axes1 = axes;
hold(axes1,'on');



% rms residual values to plot
lvls = [0, 1:5];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
[c1,h1] = contourf(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC)./rms_uncert, lvls, 'LineWidth',4,...
    'EdgeColor', 'k');
clabel(c1,h1,'FontSize',20,'FontWeight','bold');


% rms residual values to plot
% lvls = [0, 0.004, 0.005, 0.01:0.01:1];
% 
% [c1,h1] = contourf(R_bot_fine(:,:, idx_tauC), R_top_fine(:,:, idx_tauC), rms_residual(:,:, idx_tauC), lvls, 'LineWidth',4,...
%     'EdgeColor', 'k');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
    ' between EMIT and LibRadTran'],'Interpreter','latex', 'FontSize', 33);

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')


% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])



%% Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----



% define the optical depth slice you'd like to plot
idx_rTop = r_top_fine == r_top_min(1);

% Create figure
figure;


% Create axes
axes1 = axes;
hold(axes1,'on');


% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour
% [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)), n, 'LineWidth',4, 'EdgeColor', mySavedColors(9, 'fixed'));
% clabel(c1,h1,'FontSize',20,'FontWeight','bold');


% rms residual values to plot
lvls = [0, 0.003, 0.005, 0.01:0.01:1];

% Create contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :), length(r_bot_fine),...
    length(tau_c_fine)), lvls, 'LineWidth',4, 'EdgeColor', 'k');



% Create ylabel
ylabel('$r_{top}^{min} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['Relative RMS between EMIT and LibRadTran at min $r_{top}$'],'Interpreter','latex', ...
    'Fontsize', 23);

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')

% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])



%% Create a contour plot that varies with r top and r bottom and shows more than one optical depth



% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC = tau_c_fine == tau_c_min(1);
%idx_tauC = tau_c_fine == 6.2;

% Create figure
figure;


% Create axes
axes1 = axes;
hold(axes1,'on');


% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
lvls = [0, 0.3, 0.5, 1:2];


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
[c1,h1] = contour(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC)./rms_uncert, lvls, 'LineWidth',4,...
    'EdgeColor', mySavedColors(2, 'fixed'));
clabel(c1,h1,'FontSize',20,'FontWeight','bold');

% plot another optical depth
hold on
% plot the mimimum rms residual
idx_tauC = tau_c_fine == 6.3;

[c2,h2] = contour(r_bot_fine, r_top_fine, rms_residual(:,:, idx_tauC)./rms_uncert, lvls, 'LineWidth',4,...
    'EdgeColor', mySavedColors(3, 'fixed'));
clabel(c2,h2,'FontSize',20,'FontWeight','bold');

% create a legend
legend('$\tau_c = 6.4$', '$\tau_c = 6.3$', 'Interpreter', 'latex', 'Fontsize', 30', 'location', 'best')


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
    ' between EMIT and LibRadTran'],'Interpreter','latex', 'FontSize', 33);

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')


% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])


%% Create a surface plot that shows the volume of the solution space
% The set of possible solutions, the volume in (r_top, r_bot, and tau_c)
% is the space where the rms difference between the measurements and the 
% estimates are less than the rms of the measurement uncertainty


% Lets find all x, y and z values where the rms(R(x) - m)/rms(delta m) is
% less than 1


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
idx = rms_residual./rms_uncert < 0.3;

figure; 

f = fill3(R_bot_fine(idx), R_top_fine(idx), Tau_c_fine(idx), 'r');

f.EdgeAlpha = 1;

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

grid on; grid minor

% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
set(gcf, 'Position', [0 0 900 900])

xlim([r_bot_fine(1) , r_bot_fine(end)])
ylim([r_top_fine(1) , r_top_fine(end)])
zlim([tau_c_fine(1) , tau_c_fine(end)])



%% Create a surface plot at the optical depth associated with the minimum rms 
% plot the RMS residual at the minimum optical depth and let the radii at
% cloud top and bottom varry


% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC = tau_c_fine == tau_c_min(1);
%idx_tauC = tau_c_fine == 6.2;

% Create figure
figure;



% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
s = surf(R_bot_fine(:,:, idx_tauC), R_top_fine(:,:, idx_tauC), rms_residual(:,:, idx_tauC)./rms_uncert);

s.EdgeAlpha = 0.5;

% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')



% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$RMS(R(\vec{x}) - \vec{m}) / RMS(\delta \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);


% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
    ' between EMIT and LibRadTran'],'Interpreter','latex', 'FontSize', 33);

box(gca,'on');
grid(gca,'on');
axis(gca,'tight');
hold(gca,'off');
% Set the remaining axes properties
set(gca,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')


% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])




%% Make a plot of the RMS data cube showing all three variables

figure; 

scatter3(R_bot_fine(:), R_top_fine(:), Tau_c_fine(:), 40, rms_residual(:),'filled')    % draw the scatter plot(:)

% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

grid on; grid minor


set(gcf, 'Position', [0 0 1200 1200])


%%





% -----------------------------------------------------------------
% ----------- COMPARISON PLOTS BETWEEN EMIT AND MODIS -------------
% -----------------------------------------------------------------








%% Make a plot of the volume of solution space that is less than the rms uncertainty


idx_uncert = rms_residual./rms_uncert <= 1;
idx_uncert_MODIS7 = rms_residual_MODIS7./rms_uncert_MODIS7 <= 1;

percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rms_residual);

percent_states_less_than_rms_uncert_MODIS7 = sum(idx_uncert_MODIS7, 'all')/numel(rms_residual_MODIS7);


figure; 

subplot(1,2,1)

% plot volume of states using 35 EMIT measurements
s = scatter3(R_bot_fine(idx_uncert), R_top_fine(idx_uncert), Tau_c_fine(idx_uncert), 40, rms_residual(idx_uncert),'filled');    % draw the scatter plot(:)
%change face alpha
%s.MarkerFaceColor = mySavedColors(1, 'fixed');
s.MarkerFaceAlpha = 0.5;
% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')

% set limits
xlim([r_bot_fine(1), r_bot_fine(end)])
ylim([r_top_fine(1), r_top_fine(end)])
zlim([tau_c_fine(1), tau_c_fine(end)])

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['35 Channels - $RMS(R(\vec{x}) - \vec{m}) \leq RMS(\delta \vec{m})$'],'Interpreter','latex', 'FontSize', 33);

grid on; grid minor


% plot volume of states using 7 MODIS channels

subplot(1,2,2)
s2 = scatter3(R_bot_fine(idx_uncert_MODIS7), R_top_fine(idx_uncert_MODIS7), Tau_c_fine(idx_uncert_MODIS7),...
    40, rms_residual(idx_uncert_MODIS7),'filled');    % draw the scatter plot(:)

s2.MarkerFaceAlpha = 0.5;
% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')


% set limits
xlim([r_bot_fine(1), r_bot_fine(end)])
ylim([r_top_fine(1), r_top_fine(end)])
zlim([tau_c_fine(1), tau_c_fine(end)])


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['7 Channels - $RMS(R(\vec{x}) - \vec{m}) \leq RMS(\delta \vec{m})$'],'Interpreter','latex', 'FontSize', 33);

grid on; grid minor


set(gcf, 'Position', [0 0 2400 1200])




%% Create a surface plot at the optical depth associated with the minimum rms 
% plot the RMS residual at the minimum optical depth and let the radii at
% cloud top and bottom varry


% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC_emit = tau_c_fine == tau_c_min(1);
%idx_tauC = tau_c_fine == 6.2;

% Create figure
figure;



% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
s = surf(R_bot_fine(:,:, idx_tauC_emit), R_top_fine(:,:, idx_tauC_emit), rms_residual(:,:, idx_tauC_emit)./rms_uncert);

s.EdgeAlpha = 0.5;
s.FaceAlpha = 0.25;
s.FaceColor = mySavedColors(1, 'fixed');

hold on;


% --------- Plot MODIS 7 results ------------
% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC_modis = tau_c_fine == tau_c_min_MODIS7(1);

s2 = surf(R_bot_fine(:,:, idx_tauC_modis), R_top_fine(:,:, idx_tauC_modis), rms_residual_MODIS7(:,:, idx_tauC_modis)./rms_uncert_MODIS7);

s2.EdgeAlpha = 0.5;
s2.FaceAlpha = 0.25;
s2.FaceColor = mySavedColors(2, 'fixed');


% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, '$RMS(R(\vec{x}) - \vec{m}) / RMS(\delta \vec{m})$', 'FontSize', 30, 'Interpreter', 'latex')



% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$RMS(R(\vec{x}) - \vec{m}) / RMS(\delta \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);


% Create title
title(['RMS Residual between EMIT and LibRadTran'],'Interpreter','latex', 'FontSize', 33);

% define legend
legend(['35 EMIT wavelengths - min $\tau_c = $', num2str(tau_c_fine(idx_tauC_emit))],...
    ['7 EMIT wavelengths - min $\tau_c = $', num2str(tau_c_fine(idx_tauC_modis))], 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', 25)

box(gca,'on');
grid(gca,'on');
axis(gca,'tight');
hold(gca,'off');
% Set the remaining axes properties
set(gca,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')


set(gcf, 'Position', [0 0 1200 1200])




%% Code suggested by Claude

% RMS Residual Visualization Script
% Visualization of RMS residual across r_top_fine, r_bot_fine, and tau_c_fine

% Note: Replace this with your actual data loading method
% Assuming you have a 3D matrix of RMS residuals pre-calculated
% If not, you'll need to add nested loops to calculate rms_residual for each combination

% Load or generate data (MODIFY THIS SECTION)
% Example: load('your_data_file.mat', 'r_top_fine', 'r_bot_fine', 'tau_c_fine', 'rms_residual');

% If data is not pre-calculated, you might need something like:
% rms_residual = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine));
% for i = 1:length(r_top_fine)
%     for j = 1:length(r_bot_fine)
%         for k = 1:length(tau_c_fine)
%             % Calculate rms_residual for each combination
%             % rms_residual(i,j,k) = ... your calculation here
%         end
%     end
% end

% Visualization Option 1: Comprehensive 3D Visualization
figure('Position', [100, 100, 1500, 500]);

% Slice Plot
subplot(1,3,1);
slice(r_top_fine, r_bot_fine, tau_c_fine, rms_residual, ...
    mean(r_top_fine), mean(r_bot_fine), mean(tau_c_fine));
title('RMS Residual Slice Plot');
xlabel('r_{top,fine}');
ylabel('r_{bot,fine}');
zlabel('\tau_{c,fine}');
colorbar;
colormap('jet');

% Isosurface Plot
subplot(1,3,2);
p = patch(isosurface(r_top_fine, r_bot_fine, tau_c_fine, rms_residual, median(rms_residual(:))));
title('RMS Residual Isosurface');
xlabel('r_{top,fine}');
ylabel('r_{bot,fine}');
zlabel('\tau_{c,fine}');
isonormals(r_top_fine, r_bot_fine, tau_c_fine, rms_residual, p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 0.7;
lighting gouraud;
camlight;

% 2D Contour Plot (Projection at median tau_c_fine)
subplot(1,3,3);
median_tau_index = round(length(tau_c_fine)/2);
rms_2d = squeeze(rms_residual(:,:,median_tau_index));
contourf(r_top_fine, r_bot_fine, rms_2d');
title(['RMS Residual Contour (at \tau_{c,fine} = ' num2str(tau_c_fine(median_tau_index)) ')']);
xlabel('r_{top,fine}');
ylabel('r_{bot,fine}');
colorbar;
colormap('jet');

% Additional Visualization: Heatmap across different tau_c_fine values
figure;
tiledlayout(2,2);

% Create multiple heatmaps for different tau_c_fine values
tau_indices = [1, round(length(tau_c_fine)/4), round(length(tau_c_fine)/2), length(tau_c_fine)];

for i = 1:length(tau_indices)
    nexttile;
    current_rms_2d = squeeze(rms_residual(:,:,tau_indices(i)));
    heatmap(r_bot_fine, r_top_fine, current_rms_2d, ...
        'Title', ['\tau_{c,fine} = ' num2str(tau_c_fine(tau_indices(i)))], ...
        'XLabel', 'r_{bot,fine}', ...
        'YLabel', 'r_{top,fine}');
end

% Statistical Analysis
fprintf('RMS Residual Statistics:\n');
fprintf('Minimum RMS Residual: %f\n', min(rms_residual(:)));
fprintf('Maximum RMS Residual: %f\n', max(rms_residual(:)));
fprintf('Mean RMS Residual: %f\n', mean(rms_residual(:)));
fprintf('Median RMS Residual: %f\n', median(rms_residual(:)));

% Optional: Find the combination with minimum RMS residual
[min_rms, linear_index] = min(rms_residual(:));
[i, j, k] = ind2sub(size(rms_residual), linear_index);
fprintf('\nLowest RMS Residual Details:\n');
fprintf('r_{top,fine} = %f\n', r_top_fine(i));
fprintf('r_{bot,fine} = %f\n', r_bot_fine(j));
fprintf('tau_{c,fine} = %f\n', tau_c_fine(k));
fprintf('Minimum RMS Residual = %f\n', min_rms);

%%





% -------------------------------------
% ----------- MODIS PLOTS -------------
% -------------------------------------






%% Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
% ***  USING JUST FIRST 7 MODIS SPECTRAL CHANNELS ***   
% plot the RMS residual at the minimum optical depth and let the radii at
% cloud top and bottom varry


% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC = tau_c_fine == tau_c_min_MODIS7(1);

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');


% rms residual values to plot
lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];



% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
[c1,h1] = contourf(r_bot_fine, r_top_fine, rms_residual_MODIS7(:,:, idx_tauC)./rms_uncert_MODIS7, lvls, 'LineWidth',4,...
    'EdgeColor', 'k');
clabel(c1,h1,'FontSize',20,'FontWeight','bold');




% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
    ' between EMIT and LibRadTran using first 7 MODIS wavelengths'],'Interpreter','latex', 'FontSize', 20);

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');



% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])



%% Create Contour plot of rms residual between true EMIT measurements and the libRadTran modeled measurements
% --- (r_top - r_bot) versus tau  for the minimum r_top ----
% ***  USING JUST FIRST 7 MODIS SPECTRAL CHANNELS *** 


% define the optical depth slice you'd like to plot
idx_rTop = r_top_fine == r_top_min_MODIS7(1);

% Create figure
figure;


% Create axes
axes1 = axes;
hold(axes1,'on');


% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour
% [c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
%     length(tau_c_fine)), n, 'LineWidth',4, 'EdgeColor', mySavedColors(9, 'fixed'));
% clabel(c1,h1,'FontSize',20,'FontWeight','bold');


% rms residual values to plot
lvls = [0, 0.003, 0.005, 0.01:0.01:1];

% Create contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual_MODIS7(idx_rTop,:, :), length(r_bot_fine),...
    length(tau_c_fine)), lvls, 'LineWidth',4, 'EdgeColor', 'k');



% Create ylabel
ylabel('$r_{top}^{min} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title(['Relative RMS between 1st 7 MODIS Channels and LibRadTran at min $r_{top}$'],'Interpreter','latex', ...
    'Fontsize', 23);

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'tight');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')

% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])



%% Make 3D plot of contour slices along the tau dimension




% rms residual values to plot
lvls = [0, 0.25, 0.5:0.25:4];

% define the slices along the tau dimension to plot
tau_slice = [6, 6.2, 6.4, 6.6];

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'Position',[0.11068 0.11 0.7254 0.815]);
hold(axes1,'on');

s = contourslice(R_bot_fine, R_top_fine, Tau_c_fine, rms_residual_MODIS7./rms_uncert_MODIS7, [], [], tau_slice, lvls);
view(3)

hold on; plot3(r_bot_min_MODIS7(1), r_top_min_MODIS7(1), tau_c_min_MODIS7(1), 'r.', 'markersize', 20)

% Create colorbar
cb = colorbar(axes1,'Position',[0.89108 0.11 0.016 0.815]);
% create colorbar label
ylabel(cb, '$RMS(R(\vec{x}) - \vec{m})/RMS(\delta \vec{m})$', 'FontSize', 30, 'Interpreter', 'latex')

% set all contour lines to a certain thickness
for nn = 1:length(s)
    s(nn).LineWidth = 4;
end


% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
title('Using 7 MODIS Spectral Measurements','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Set position and size of figure
set(gcf, 'Position', [0 0 1000 1000])

% set 3D point of view and grid lines
view(axes1,[-56.6186981596273 23.687116207586]);
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on');



%% Create a surface plot at the optical depth associated with the minimum rms 
% plot the RMS residual at the minimum optical depth and let the radii at
% cloud top and bottom varry


% define the optical depth slice you'd like to plot
% plot the mimimum rms residual
idx_tauC = tau_c_fine == tau_c_min_MODIS7(1);
%idx_tauC = tau_c_fine == 6.2;

% Create figure
figure;


% compute the rms of the EMIT reflectance uncertainty
Refl_emit_uncertainty = 0.05 .* Refl_emit;
rms_uncert = sqrt(mean(Refl_emit_uncertainty.^2));


% rms residual values to plot
%lvls = [0, 0.25, 0.5, 1:3];
%lvls = [0, 0.3, 0.5, 1:2];


% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
s = surf(R_bot_fine(:,:, idx_tauC), R_top_fine(:,:, idx_tauC), rms_residual_MODIS7(:,:, idx_tauC));

s.EdgeAlpha = 0.5;

% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')



% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
zlabel('$RMS(R(\vec{x}) - \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);


% Create title
title(['RMS Residual for $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
    ' between first 7 MODIS channels and LibRadTran'],'Interpreter','latex', 'FontSize', 33);

box(gca,'on');
grid(gca,'on');
axis(gca,'tight');
hold(gca,'off');
% Set the remaining axes properties
set(gca,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
    'on');
% % Create colorbar
% cb = colorbar(axes1);
% % create colorbar label
% ylabel(cb, '$1/sr$', 'FontSize', 30, 'Interpreter', 'latex')


% set the figure size to be proportional to the length of the r_top and
% r_bot vectors
%set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
set(gcf, 'Position', [0 0 900 900])

