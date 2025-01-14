

clear variables


%% LOAD DATA SET

% load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-26-Nov-2024_rev1.mat')
%load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-27-Nov-2024_rev1.mat')
load('reflectance_calcs_EMIT-data-from-17_Jan_2024_coast_sim-ran-on-11-Dec-2024_rev1.mat')


%% Define synthetic model data


r_top_truth = 10.17;
r_bot_truth = 4.74;
tau_c_truth = 5.96;

[~, idx_r_top] = min(abs(r_top_fine - r_top_truth));
[~, idx_r_bot] = min(abs(r_bot_fine - r_bot_truth));
[~, idx_tau_c] = min(abs(tau_c_fine - tau_c_truth));

synthetic_measurement = reshape(Refl_model_fine(idx_r_top, idx_r_bot, idx_tau_c, :), [],1);


%% Add Gaussian Noise to the measurements

% --- meausrement uncertainty ---
% define this as a fraction of the measurement
measurement_uncert = 0.01;

% Define a gaussian where the mean value is the true measurement, and twice
% the standard deviation is the product of the measurement uncertainty and
% the true measurements. 
% Remember: +/- 1*sigma = 68% of the area under the gaussian curve
%           +/- 2*sigma = 95% of the area under the gaussian curve
%           +/- 3*sigma = 99.7% of the area under the gaussian curve

% Compute the new synethtic measurement with gaussian noise
% *** Gaussian noise can be either positive or negative. Meaning, an
% uncertainty of 5% implies the true value can lie anywhere between the
% measured value +/- 5% of the measured value
% define a 
synthetic_measurement_with_noise = synthetic_measurement + synthetic_measurement.*(measurement_uncert/3) .*...
    randn(length(inputs.bands2run), 1);

% define the synthetic relfectance uncertainty
synthetic_measurement_uncert = measurement_uncert .* synthetic_measurement_with_noise;




%%



% % Create the new fine grid to interpolate on
% % define the discrete step length of each variable
% d_r_top = 0.1;      % microns
% d_r_bot = 0.1;      % microns
% d_tau_c = 0.1;
% 
% r_top_fine = r_top(1):d_r_top:r_top(end);
% r_bot_fine = r_bot(1):d_r_bot:r_bot(end);
% tau_c_fine = tau_c(1):d_tau_c:tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

% Let's now seperate out the interpolated relfectance at the seven MODIS
% wavelengths
wl_MODIS7_idx = [1, 4, 6, 7, 19, 23, 29];

% define the synthetic measurement
synthetic_measurement_with_noise_MODIS7 = synthetic_measurement_with_noise(wl_MODIS7_idx);
synthetic_measurement_uncert_MODIS7 = synthetic_measurement_uncert(wl_MODIS7_idx);

% Grab the modeled data at just the 7 MODIS bands
Refl_model_fine_MODIS7 = Refl_model_fine(:,:,:, wl_MODIS7_idx);




%% Using an exact modeled estimate without noise

use_l2_norm = true;


if use_l2_norm==false

    % Compute the rms difference between the measurements and the modeled
    % reflectances
    rms_residual = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    rms_residual_MODIS7 = sqrt( mean( (repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine_MODIS7).^2, 4));

    % Compute the RMS of the synthetic measurement uncertainty
    rms_uncert = sqrt( mean( synthetic_measurement_uncert.^2));

    % Compute the RMS of the synthetic measurement uncertainty using 7
    % MODIS channels
    rms_uncert_MODIS7 = sqrt( mean( synthetic_measurement_uncert_MODIS7.^2));


else

    %Compute the l2 norm difference between the measurements and the modeled reflectances
    rms_residual = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine).^2, 4));

    rms_residual_MODIS7 = sqrt( sum( (repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine))...
        - Refl_model_fine_MODIS7).^2, 4));

    % Compute the l2 norm of the synthetic measurement uncertainty
    rms_uncert = sqrt( sum( synthetic_measurement_uncert.^2));

    % Compute the l2 norm of the synthetic measurement uncertainty using 7
    % MODIS channels
    rms_uncert_MODIS7 = sqrt( sum( synthetic_measurement_uncert_MODIS7.^2));


end


% Using the new fine grid, calculate how many sets of measurements are
% within the EMIT measurement and it's uncertainty
redundant_states = all( abs( repmat(reshape(synthetic_measurement_with_noise, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine)) -...
    Refl_model_fine) <= repmat(reshape(synthetic_measurement_uncert, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine)) ,4);

redundant_states_MODIS7 = all( abs( repmat(reshape(synthetic_measurement_with_noise_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine)) -...
    Refl_model_fine_MODIS7) <= repmat(reshape(synthetic_measurement_uncert_MODIS7, 1, 1, 1, []), length(r_top_fine), length(r_bot_fine), length(tau_c_fine)) ,4);


%% Find the states with the lowest rms residul

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
s = surf(R_bot_fine(:,:, idx_tauC), R_top_fine(:,:, idx_tauC), rms_residual(:,:, idx_tauC));

s.EdgeAlpha = 0.5;

% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')


% -----------------------------------------------------------------------
% show the fraction of points where the rms residual is less than the rms
% of the measurement uncertainty
% idx_solution_space = rms_residual <= rms_uncert;
% hold on; plot3(R_bot_fine(idx_solution_space), R_top_fine(idx_solution_space), rms_residual(idx_solution_space),...
%     '.k', 'MarkerSize', 20)

% try ploting a transparent plane
hold on; surf(R_bot_fine(:,:, idx_tauC), R_top_fine(:,:, idx_tauC), repmat(rms_uncert, length(r_top_fine), length(r_bot_fine)),...
    'FaceColor', 'k', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.25)
% -----------------------------------------------------------------------

% Create ylabel
ylabel('$r_{top}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);




% Create title
if use_l2_norm==false

    % Create xlabel
    zlabel('$RMS(R(\vec{x}) - \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    title(['RMS Residual at global min $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

else

    % Create xlabel
    zlabel('$RSS(R(\vec{x}) - \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    title(['RSS Residual at global min $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

end



zlim([0, max(rms_residual, [], 'all')])

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
set(gcf, 'Position', [0 0 1150 1100])


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
if use_l2_norm==false

    title(['RMS Residual at global min $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

else

    title(['RSS Residual at global min $\tau_c = $', num2str(tau_c_fine(idx_tauC)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

end



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
lvls = [0, 1:24];


% Create contour
[c1,h1] = contourf(tau_c_fine, r_top_min(1)-r_bot_fine, reshape(rms_residual(idx_rTop,:, :)./rms_uncert, length(r_bot_fine),...
    length(tau_c_fine)),  lvls, 'LineWidth',4, 'EdgeColor', 'k');
clabel(c1,h1,'FontSize',20,'FontWeight','bold');




% Create ylabel
ylabel('$r_{top}^{min} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create title
if use_l2_norm==false

    title(['RMS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

else

    title(['RSS Residual at global min $r_{top} = $', num2str(r_top_fine(idx_rTop)),...
        ' between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

end


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



%% Create a surface plot where the x-axis is (r_top - b_bot), the y-axis is optical depth and the 
% z-axis is the rms residual




% compute the difference between r_top and r_bot at each grid point
R_top_minus_R_bot = R_top_fine(:,:,1) - R_bot_fine(:,:,1);

% stack columns
R_top_minus_R_bot = reshape(R_top_minus_R_bot, [], 1);

% Let's sort it in ascending order
[R_top_minus_R_bot_sort, idx_sort] = sort(R_top_minus_R_bot);


% Let's collapse the rms_residual matrix into two dimentions
% we need to reshape the rms_residaul matrix

% r_top changes along the y-dimension, or the row-space
% r_bot changes along the x-dimension, or the column-space
% tau_c changes along the z-dimension, or the depth-space

% --- Stack columns just like above! ----
%rms_residual_2D = reshape(rms_residual, [], size(rms_residual, 3));

rms_residual_2D = zeros(size(rms_residual, 1)*size(rms_residual, 2), size(rms_residual, 3));

for dd = 1:size(rms_residual, 3)

    data_column = [];

    for cc = 1:size(rms_residual, 2)
    
        data_column = [data_column; rms_residual(:,cc,dd)];

    end

    % Sort the data column in ascending order of the independent variable
    rms_residual_2D(:,dd) = data_column(idx_sort);

end






% Create a mesh grid!
R_top_minus_R_bot_sort = repmat(R_top_minus_R_bot_sort, 1, size(rms_residual, 3));

% create a matrix with the same dimensions as R_bot_minus_R_top but with
% optical depth changing along row spac
Tau_c_2D = repmat(tau_c_fine, size(R_top_minus_R_bot,1), 1);




% Create a surface plot
figure;

% % Create contour plot showing all radii at cloud top and bottom for a
% % particular optical depth
% s = surf(R_top_minus_R_bot_sort, Tau_c_2D, rms_residual_2D);


% only plot the states where the rms residual is less than the rms of the
% measurement uncertainty

idx_uncert = rms_residual_2D <= rms_uncert;

% Create contour plot showing all radii at cloud top and bottom for a
% particular optical depth
s = surf(R_top_minus_R_bot_sort(idx_uncert), Tau_c_2D(idx_uncert), rms_residual_2D(idx_uncert));


s.EdgeAlpha = 0.5;

% Create colorbar
cb = colorbar;
% create colorbar label
ylabel(cb, 'Reflectance ($1/sr$)', 'FontSize', 30, 'Interpreter', 'latex')


% Create ylabel
ylabel('$\tau_c$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

% Create xlabel
xlabel('$r_{top} - r_{bot}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);




% Create title
if use_l2_norm==false

    % Create xlabel
    zlabel('$RMS(R(\vec{x}) - \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    title(['RMS Residual between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

else

    % Create xlabel
    zlabel('$RSS(R(\vec{x}) - \vec{m})$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    title(['RSS Residual between Synthetic Measurements with ', num2str(100*measurement_uncert),...
        '\% uncertainty and LibRadTran'],'Interpreter','latex', 'FontSize', 23);

end



zlim([0, max(rms_residual, [], 'all')])

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
set(gcf, 'Position', [0 0 1150 1100])



%% 

idx_uncert = rms_residual./rms_uncert <= 1;
idx_uncert_MODIS7 = rms_residual_MODIS7./rms_uncert_MODIS7 <= 1;

percent_states_less_than_rms_uncert = sum(idx_uncert, 'all')/numel(rms_residual);

percent_states_less_than_rms_uncert_MODIS7 = sum(idx_uncert_MODIS7, 'all')/numel(rms_residual_MODIS7);





%% plots from Claude

% RMS Residual Visualization Script
% Visualization of RMS residual across r_top_fine, r_bot_fine, and tau_c_fine

% Note: Replace this with your actual data loading method
% Assuming you have a 3D matrix of RMS residuals pre-calculated
% If not, you'll need to add nested loops to calculate rms_residual for each combination

% Load or generate data (MODIFY THIS SECTION)
% Example: load('your_data_file.mat', 'r_top_fine', 'r_bot_fine', 'tau_c_fine', 'rms_residual');

% Visualization Option 1: Comprehensive 3D Visualization
figure('Position', [100, 100, 1500, 1500]);

% Slice Plot
subplot(1,3,1);
slice(r_bot_fine, r_top_fine, tau_c_fine, rms_residual, ...
    mean(r_bot_fine), mean(r_top_fine), mean(tau_c_fine));
title('RMS Residual Slice Plot');
xlabel('r_{bot}');
ylabel('r_{top}');
zlabel('\tau_{c}');
colorbar;
colormap('jet');

% Isosurface Plot
subplot(1,3,2);
p = patch(isosurface(r_bot_fine, r_top_fine, tau_c_fine, rms_residual, median(rms_residual(:))));
title('RMS Residual Isosurface');
xlabel('r_{bot}');
ylabel('r_{top}');
zlabel('\tau_{c}');
isonormals(r_bot_fine, r_top_fine, tau_c_fine, rms_residual, p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 0.7;
lighting gouraud;
camlight;

% 2D Contour Plot (Projection at median tau_c_fine)
subplot(1,3,3);
median_tau_index = round(length(tau_c_fine)/2);
rms_2d = squeeze(rms_residual(:,:,median_tau_index));
contourf(r_bot_fine, r_top_fine, rms_2d);
title(['RMS Residual Contour (at \tau_{c,fine} = ' num2str(tau_c_fine(median_tau_index)) ')']);
xlabel('r_{bot}');
ylabel('r_{top}');
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
