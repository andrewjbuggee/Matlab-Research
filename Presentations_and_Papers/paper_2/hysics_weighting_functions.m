%% Weighting functions from the 66 HySICS channels

clear variables

% 66 HySICS bands
load('disort_HySICS_reflectance_for_weightingFunctions_inhomogeneous_droplet_profile_sim-ran-on-17-Jul-2025_rev1.mat')


% 7 MODIS bands
% load('disort_HySICS_reflectance_for_weightingFunctions_inhomogeneous_droplet_profile_MODIS-7bands_sim-ran-on-17-Jul-2025_rev1.mat')

% ********************************************
% *** Vary The Optical Thickness Linearly! ***
% ********************************************
% tau_2run = linspace(0.001, inputs.RT.tau_c, size(Refl_model,1))';


% length of each independent variable
num_wl = length(inputs.bands2run);
num_tau_layers = length(tau_2run);

num_INP_files = num_wl*num_tau_layers;


% Let's fit a moving average to each weighting function and the renormalize

f = zeros(size(w));

N_mov_avg = 30;


tau_midPoint = tau_2run(1:end-1,:) + diff(tau_2run, 1, 1);

for ww = 1:num_wl

    % find the moving average
    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    %f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);
    f(:,ww) = movmean(w(:,ww), N_mov_avg);


    % renormalize!
    a = 1/trapz(tau_midPoint, f(:,ww));

    f(:,ww) = f(:,ww).*a;


end


% plot just the moving average weighting functions


% *** define which wavelengths to plot ***
% Full set of all 66 HySICS channels:
% wl_2plot = inputs.RT.wavelengths2run(:,1);
%
% Paper 2, Figure 1 uses the "best" 17-wavelength subset of the 66 (one
% visible channel plus the 940 nm / 1.9-2.1 um water-vapor absorption
% complex). These are all exact members of the 66-channel set.
wl_2plot = [540 1432 1490 1766 1934 1937 1940 1943 1952 1956 1962 ...
    1968 1974 1980 1992 2005 2075]';

lgnd_str = cell(numel(wl_2plot), 1);

% *** colorblind-friendly line colors ***
% cividis is a perceptually-uniform sequential colormap designed for
% color-vision deficiency (Nunez, Anderton & Renslow, 2018). One color per
% line, ordered by wavelength, so the color encodes the wavelength ordering.
% Generate 2 extra colors and drop the palest so the brightest lines stay
% legible on a white background.
% --- navy to yellow map ---
% cmap = cividis(numel(wl_2plot) + 2);
% cmap = cmap(1:end-2, :);

% --- navy to pink to yellow to white map ---
cmap = ametrine(numel(wl_2plot) + 2);
cmap = cmap(1:end-2, :);

% cmap = mySavedColors(61:(61+numel(wl_2plot)-1), 'fixed');

fig1 = figure;


if inputs.RT.monochromatic_calc==true

    for ww = 1:length(wl_2plot)

        [~,idx2plot] = min(abs(inputs.RT.wavelengths2run(:,1) - wl_2plot(ww)));

        plot(f(:,idx2plot), tau_midPoint, 'LineStyle', '-', ...
            'Color', cmap(ww,:), 'LineWidth', 2.5)

        hold on

        lgnd_str{ww} = ['$\lambda = $', num2str(round(inputs.RT.wavelengths2run(idx2plot, 1))), ' $nm$'];

    end



else

    tau = changing_variables(2:end,3);

    % plot(w, diff(flipud(changing_variables(:,2)))/2 + flipud(changing_variables(2:end,2)))
    plot(w, flipud(tau))

    % --- overlay a smoothed spline fit ---
    % Create smooth spline function
    % f=fit(diff(flipud(changing_variables(:,2)))/2 + flipud(tau), w, 'smoothingspline','SmoothingParam',0.95);
    f = movmean(w, N_mov_avg);

end




% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$w(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
% title(['Weighting Function at ', num2str(changing_variables(1,1)), ' nm'],'Interpreter','latex')
% title('Weighting Functions for first 7 MODIS Spectral Channels','Interpreter','latex')




set(gcf, 'Position',[0 0 1400 800])




% Create Legend
% legend(string(inputs.RT.wavelengths2run(:,1))','Interpreter','latex','Location','northwest','FontSize',22)
legend(lgnd_str,'Interpreter','latex','Location','eastoutside','FontSize',22)

% Create textbox with simulation properties

% Textbox
dim = [0.155714285714286 0.144548492431641 0.196462309701102 0.382951507568359];

if ischar(inputs.RT.sensor_altitude)==true
    texBox_str = {['$sza$ = ',num2str(inputs.RT.sza)],...
        ['$vza$ = ',num2str(inputs.RT.vza)],...
        ['$z_{out}$ = ', inputs.RT.sensor_altitude],...
        ['$Cloud\;top$ = ', num2str(inputs.RT.z_topBottom(1)), ' km'],...
        ['$Cloud\;base$ = ', num2str(inputs.RT.z_topBottom(2)), ' km'],...
        ['$r_{top}$ = ',num2str(round(inputs.RT.r_top)), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.RT.r_bot)), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.RT.tau_c)],...
        ['$A_0$ = ', num2str(inputs.RT.surface_albedo)]};
else
    texBox_str = {['$sza$ = ',num2str(inputs.RT.sza)],...
        ['$vza$ = ',num2str(inputs.RT.vza)],...
        ['$z_{out}$ = ', num2str(inputs.RT.sensor_altitude), ' km'],...
        ['$Cloud\;top$ = ', num2str(inputs.RT.z_topBottom(1)), ' km'],...
        ['$Cloud\;base$ = ', num2str(inputs.RT.z_topBottom(2)), ' km'],...
        ['$r_{top}$ = ',num2str(round(inputs.RT.r_top)), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.RT.r_bot)), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.RT.tau_c)],...
        ['$A_0$ = ', num2str(inputs.RT.surface_albedo)]};
end

t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig1,[folderpath_figs,'HySICS Weighting Functions.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(fig1,[folderpath_figs,'HySICS weighting functions.jpg'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------
