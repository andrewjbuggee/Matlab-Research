%% Create Plots for IRS 2022 Presentation!


% By Andrew John Buggee


%% Create a plot for of the single scattering albedo of a water droplet

justQ = false;
dist_str = 'mono';              % mono dispersed distribution
wl = 100:2300;                  % wavelengths (nm)
re = [1, 10:10:100];            % microns - effective radius

yq = zeros(length(wl), 8, length(re));

f = figure;

for rr = 1:length(re)
    yq(:,:,rr) = interp_mie_computed_tables([wl',linspace(re(rr), re(rr),length(wl))'],dist_str,justQ);

    % Plot the single scattering albedo
    omega = yq(:,6,rr);
    plot(wl,omega(:),'LineWidth',3); hold on
end

xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Single Scattering Albedo', 'Interpreter','latex')
grid on; grid minor
legend(strcat(string(re), ' $\mu m$'),'Location','best','Interpreter','latex')


%%