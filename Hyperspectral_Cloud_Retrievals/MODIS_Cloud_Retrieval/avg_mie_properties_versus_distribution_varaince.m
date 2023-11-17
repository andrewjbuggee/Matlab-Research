%% How do single scattering properties averaged over a size distribution change with distribution variance?

clear variables

% By Andrew John Buggee

%% Define the cloud parameters that will be changing during each reflectance calculation

% define the type of droplet distribution
distribution_str = 'gamma';

% define the distribution varaince
dist_var = [1, 100];

% define the droplet effective radius
re = 5;                      % microns

% set up a wavelength range to compute mie properties
wavelength_min = 350;           % nm
wavelength_max = 2500;           % nm
wavelength_step = 2;            % nm

index_of_refraction = 'water';


%% compute the average properties

% Define the wavelength grid
wavelength = [wavelength_min, wavelength_max, wavelength_step];
% define the wavelength vector
wl_vec = (wavelength_min:wavelength_step:wavelength_max)';


% set up the empty arrays
ssa_avg = zeros(numel(wl_vec), numel(dist_var));
Qe_avg = zeros(numel(wl_vec), numel(dist_var));
g_avg = zeros(numel(wl_vec), numel(dist_var));
legend_str = cell(1, numel(dist_var));


for vv = 1:numel(dist_var)

    [ssa_avg(:,vv), Qe_avg(:,vv), g_avg(:,vv)] = average_mie_over_size_distribution(re, dist_var(vv), wavelength,...
        index_of_refraction, distribution_str);

    legend_str{vv} = ['$r_e = $', num2str(re), '$\mu m$, $\sigma = $', num2str(dist_var(vv))];


end


%% Plot the resutls

markerSize = 20;


% --- first plot the average single scattering albedo ---
figure1 = figure('WindowState','maximized');
% Create axes
axes1 = axes('Parent',figure1,'Position',[0.08125 0.31578 0.25 0.5]);
hold(axes1,'on');
plot(wl_vec, ssa_avg, '.', 'MarkerSize', markerSize)

grid on; grid minor

xlabel('$\lambda$ ($nm$)','Interpreter','latex', 'FontSize',35)
ylabel('$\varpi_0$','Interpreter','latex', 'FontSize',35)
axis square
box on


% --- Second plot the average extinction efficiency ---
axes2 = axes('Parent',figure1,'Position',[0.40390 0.31578 0.25 0.5]);
hold(axes2,'on');
plot(wl_vec, Qe_avg, '.', 'MarkerSize', markerSize)

grid on; grid minor

xlabel('$\lambda$ ($nm$)','Interpreter','latex', 'FontSize',35)
ylabel('$Q_e$','Interpreter','latex', 'FontSize',35)
axis square
box on


% --- Last plot the average asymmetry parameter ---
axes3 = axes('Parent',figure1,'Position',[0.71914 0.31578 0.25 0.5]);
hold(axes3,'on');
plot(wl_vec, g_avg, '.', 'MarkerSize', markerSize)

grid on; grid minor

xlabel('$\lambda$ ($nm$)','Interpreter','latex', 'FontSize',35)
ylabel('$g$','Interpreter','latex', 'FontSize',35)
axis square
box on


% write the legend, where each curve has a different effective raidus
legend(legend_str, 'Interpreter', 'latex', 'location','best', 'Fontsize',20)

% set the figure size
set(gcf, 'Position', [0 0 1000 600])




%% What do the averaged values in the LibRadTran pre-computed mie table look like?










