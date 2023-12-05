%% 3D Interpolate the reflectance calculations on a finer grid and compare with MODIS measurement
% This will plot the state space that results in an RMS residual less than
% the MODIS measurement uncertainty RMS 

function [] = plot_state_space_lessThan_MODIS_measurement_uncertainty(r_bot, r_top, tau_c, R_model, modis,...
    idx_modis_pixel)

% define the number of spectral channels used
num_bands = size(modis.EV1km.reflectance,3);

% Meshgrid is defined on x,y,z space, not row, column, depth space
% In 3D space, z = row, x = column, y = depth
[R_bot, R_top, Tau_c] = meshgrid(r_bot, r_top, tau_c);

% Create the new fine grid to interpolate on
% define the discrete step length of each variable
d_r_top = 0.1;      % microns
d_r_bot = 0.1;      % microns
d_tau_c = 0.1;

r_top_fine = r_top(1):d_r_top:r_top(end);
r_bot_fine = r_bot(1):d_r_bot:r_bot(end);
tau_c_fine = tau_c(1):d_tau_c:tau_c(end);

[R_bot_fine, R_top_fine, Tau_c_fine] = meshgrid(r_bot_fine, r_top_fine, tau_c_fine);

R_model_fine = zeros(length(r_top_fine), length(r_bot_fine), length(tau_c_fine), size(R_model,4));


% ***** IF USING REFLECTANCE CALCS USING STANDARD REFLECTANCE DEFINITION *****
% MULTIPLY EACH VALUE BY THE AIRMASS: COS(SZA)
warning([newline, 'Are you using reflectance_calcs_standardReflectance_with_mu0_9-nov-2008-data_15-Nov-2023.mat?',...
    newline, 'Make sure you multiply all reflectances by cos(sza)!', newline])


for wl = 1:size(R_model,4)

    R_model_fine(:,:,:,wl) = interp3(R_bot, R_top, Tau_c, R_model(:, :, :, wl),...
        R_bot_fine, R_top_fine, Tau_c_fine);


end



% Using the new fine grid, calculate how many sets of measurements are
% within the MODIS measurement and it's uncertainty

% Grab the MODIS reflectances for the pixel used
[r,c] = ind2sub(size(modis.EV1km.reflectance(:,:,1)), idx_modis_pixel);
R_modis = zeros(1, num_bands);
R_uncert_modis = zeros(1, num_bands);

for bb = 1:num_bands

    % ****** DID YOU USE REFLECTANCE_4MODIS? ******
    % If not you need to divide the MODIS reflectance by cos(sza)
    R_modis(bb) = modis.EV1km.reflectance(r,c,bb);
    R_uncert_modis(bb) = R_modis(bb) * 0.01*double(modis.EV1km.reflectanceUncert(r,c,bb)); % converted from percentage to reflectance
end

modis_rms_uncertainty = sqrt( 1/num_bands * sum(R_uncert_modis.^2));


states_lessThan_modisUncert = zeros(size(R_model_fine,1), size(R_model_fine,2), size(R_model_fine,3));
rms_residual = zeros(size(R_model_fine,1), size(R_model_fine,2), size(R_model_fine,3));

for rt = 1:size(R_model_fine,1)


    for rb = 1:size(R_model_fine,2)


        parfor tc = 1:size(R_model_fine,3)

            % Compute the RMS residual between the estimated measurement
            % and the true MODIS measurement
            rms_residual(rt, rb, tc) = sqrt( 1/num_bands * sum((R_modis - reshape(R_model_fine(rt,rb,tc,:), 1, [])).^2));

            % Check to see if the rms_residual is less than the MODIS
            % measurement uncertainty rms
            states_lessThan_modisUncert(rt,rb,tc) = rms_residual(rt, rb, tc) < modis_rms_uncertainty;

            


        end
    end
end

% print the number of states that have an RMS residual less than the MODIS
% uncertainty residual
if sum(states_lessThan_modisUncert, 'all')==0
    
    % if there are no measurements within the MODIS measurement
    % uncertainty, find the 500 smallest rms residuals and the associated
    % state vectors
    N = 500;
    smallest_rms_residual = zeros(1, N);
    idx_smallest = zeros(1, N);

    for nn = 1:N
        [smallest_rms_residual(nn), idx_smallest(nn)] = min(rms_residual, [], "all");
        % after grabing the smalleset value, set the value to nan
        rms_residual(idx_smallest(nn)) = nan;

    end

    [r_2plot, c_2plot, d_2plot] = ind2sub(size(rms_residual), idx_smallest);


    disp([newline, 'There are no state vectors that lead to a set of measurements that ',...
        'have a residual less than the MODIS uncertainty at each wavelength!', newline,newline,...
        'The state vector that leads to a set of measurements closest to the ',...
        'the MODIS measurements is:', newline, 'r_top = ',num2str(R_top_fine(idx_smallest(1))), '  r_bot = ',...
        num2str(R_bot_fine(idx_smallest(1))), '  tau_c = ', num2str(Tau_c_fine(idx_smallest(1))), newline,...
        'This state vector led to an rms residual of ', num2str(smallest_rms_residual(1)), newline,...
        'Plotting the state space with the smallest RMS rsidual...', newline])

    title_str = 'State space with smallest RMS residual';
    legend_str = 'Region with smallest RMS residual';

else
    disp([newline, 'There are ', num2str(sum(states_lessThan_modisUncert, 'all')), ' sets of measurements within',...
        ' the MODIS measurement and uncertainty.', newline])

    [r_2plot, c_2plot, d_2plot] = ind2sub(size(states_lessThan_modisUncert), find(states_lessThan_modisUncert));
    
    title_str = 'State space with rms residuals less than the MODIS rms uncertainty';
    legend_str = 'State Space region with RMS $<$ MODIS uncertinaty RMS';

end


% -------------------------------------------------
% ----- Plot all the states on a scatter plot -----
% -------------------------------------------------

% Lets define the color of each marker to be associated with the droplet
% size
% set the number of colors to be the length of the data to plot
r_top_2plot = zeros(length(r_2plot), 1);
r_bot_2plot = zeros(length(r_2plot), 1);
tau_c_2plot = zeros(length(r_2plot), 1);

for nn = 1:length(r_2plot)
    r_top_2plot(nn) = R_top_fine(r_2plot(nn), c_2plot(nn), d_2plot(nn));
    r_bot_2plot(nn) = R_bot_fine(r_2plot(nn), c_2plot(nn), d_2plot(nn));
    tau_c_2plot(nn) = Tau_c_fine(r_2plot(nn), c_2plot(nn), d_2plot(nn));

end

C = colormap(parula(length(tau_c_2plot)));
% sort the droplet size values
[tau_c_2plot_sort, idx_sort] = sort(tau_c_2plot, 'ascend');

figure;

for nn = 1:length(tau_c_2plot_sort)

    plot(r_bot_2plot(idx_sort(nn)), r_top_2plot(idx_sort(nn)),'Marker','.','Color',C(nn,:),'MarkerSize',25)

    hold on

end

% Plot a one-to-one line to show the boundary for homogenous profiles
[min_radius_val, ~] = min([r_top_2plot; r_bot_2plot]);
[max_radius_val, ~] = max([r_top_2plot; r_bot_2plot]);
plot([min_radius_val, max_radius_val], [min_radius_val, max_radius_val], 'k-', ...
    'linewidth', 1)

xlim([min(r_bot_2plot), max(r_bot_2plot)])
ylim([min(r_top_2plot), max(r_top_2plot)])

% set the colorbar limits
% set the limits of the colormap to be the min and max value
cb = colorbar;
if min(tau_c_2plot_sort) == max(tau_c_2plot_sort)
    clim([0.95, 1.05] * min(tau_c_2plot_sort));
else
    clim([min(tau_c_2plot_sort), max(tau_c_2plot_sort)])
end
% set colorbar title
cb.Label.String = '$\tau_c$ ($\mu m$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

grid on; grid minor
xlabel('$r_{bot}$ $(\mu m)$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
set(gcf, 'Position', [0 0 1000 500])

% Set title as the resolution of each variable
title(['$\triangle r_{top} = $', num2str(d_r_top), ' $\mu m$',...
    '    $\triangle r_{bot} = $', num2str(d_r_bot), ' $\mu m$',...
    '    $\triangle \tau_{c} = $', num2str(d_tau_c)], ...
    'Fontsize', 25, 'Interpreter', 'latex')


% ------------------------------------------------------------------------
% ---- plot the 2D space or r-top and r-bot showing area of redundancy ---
% ------------------------------------------------------------------------

% for every r-top, what is the largest and smallest r_bot that results in a
% measurement within the MODIS measurement and uncertainty?
[r_top_unique, idx_original] = unique(r_top_2plot);

top_boundary = zeros(length(r_top_unique), 2);
bottom_boundary = zeros(length(r_top_unique), 2);
tau_c_points_top = cell(length(r_top_unique),1);
tau_c_points_top_minVal = zeros(length(r_top_unique), 1);
tau_c_points_bottom = cell(length(r_top_unique),1);
tau_c_points_bottom_minVal = zeros(length(r_top_unique), 1);

for nn = 1:length(r_top_unique)

    % find set of r_bottom values for each unique r_top
    r_bottom_set = r_bot_2plot(r_top_2plot==r_top_unique(nn));

    % If there is more than 1 value in the set, find the highest and lowest
    % value. These make up the upper and lower boundaries, respectively
    % The locations should be (r_bot,r_top)
    bottom_boundary(nn,:) = [min(r_bottom_set), r_top_unique(nn)];
    top_boundary(nn,:) = [max(r_bottom_set), r_top_unique(nn)];

    % grab the optical depth for each of these points
    tau_c_points_bottom{nn} = tau_c_2plot(r_top_2plot==r_top_unique(nn) & r_bot_2plot==min(r_bottom_set));
    tau_c_points_bottom_minVal(nn) = min(tau_c_points_bottom{nn});

    tau_c_points_top{nn} = tau_c_2plot(r_top_2plot==r_top_unique(nn) & r_bot_2plot==max(r_bottom_set));
    tau_c_points_top_minVal(nn) = min(tau_c_points_top{nn});

end

% Create a polyshape using the vertices above
figure;
p = patch([bottom_boundary(:,1); flipud(top_boundary(:,1))], ...
    [bottom_boundary(:,2); flipud(top_boundary(:,2))],...
    [tau_c_points_bottom_minVal; tau_c_points_top_minVal], 'FaceColor','interp');


xlim([r_bot(1), r_bot(end)])
ylim([r_top(1), r_top(end)])

% create a one-to-one line to delineate between profiles where r-top>r-bot
% and those where this isn't true
hold on;
plot([r_top(1), r_bot(end)], [r_top(1), r_bot(end)], 'k-', 'linewidth', 2)

p.EdgeAlpha = 0;

cb = colorbar;
% set colorbar title
cb.Label.String = '$\tau_c$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 25;

% create legend
legend(legend_str, 'Vertically homogenous droplet profile',...
    'Interpreter', 'latex', 'Fontsize', 18, 'Location', 'best')

title(title_str, ...
    'Fontsize', 23)

grid on; grid minor

xlabel('$r_{bot}$ $(\mu m)$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
set(gcf, 'Position', [0 0 1200 600])




% ------------------------------------------------------------------------
% ---- Plot contours of rms residual over the state space ----------------
% ------------------------------------------------------------------------
idx_bot  = 4;

figure;

contour(reshape(Tau_c_fine(:, idx_bot, :), length(r_top_fine), length(Tau_c_fine)),...
    reshape(R_top_fine(:, idx_bot, :), length(r_top_fine), length(Tau_c_fine)),...
    reshape(rms_residual(:, idx_bot, :), length(r_top_fine), length(Tau_c_fine)), ...
    "ShowText", true, "LineWidth", 3)


grid on; grid minor

title('RMS Residual')
xlabel('$\tau_c$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
set(gcf, 'Position', [0 0 1200 600])

r_bot_slice = [];
r_top_slice = [];
tau_c_slice = [tau_c_2plot_sort(1), tau_c_2plot_sort(round(length(tau_c_2plot_sort)/2)), tau_c_2plot_sort(end)];

figure;

c = contourslice(R_bot_fine, R_top_fine, Tau_c_fine, rms_residual, r_bot_slice, r_top_slice,...
    tau_c_slice);
view(3)
set(c, 'Linewidth', 3)


grid on; grid minor

title('RMS Residual')
xlabel('$r_{bot}$ $(\mu m)$','Interpreter','latex')
ylabel('$r_{top}$ $(\mu m)$','Interpreter','latex')
zlabel('$\tau_c$','Interpreter','latex')
set(gcf, 'Position', [0 0 1200 600])

colormap hot



end
