%% Plots for my Email to Berhard Mayer


% Andrew John Bugee

%%  When using LibRadTran's Mie solver, if one assumes a droplet distribution, the calcualtions of extinction efficiency
% tend towards infinity as the effective radius tends towards 0

clear variables

%% ----- Run Mie File -----


folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_B-Mayer_mono_dispersed.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds,headers,num_radii] = readMIE(folderName,inputName(1:end-4));




%%  Plot the size parameter across all radii and wavelength values

figure; 
imagesc(ds.r_eff, ds.wavelength, 2*pi*ds.r_eff./(ds.wavelength./1e3));
colormap hot
% plot the extinction efficiency
colorbar
ylabel('Wavelength (nm)', 'Interpreter','latex')
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
title('Size Parameter - $x = \frac{2 \pi r}{\lambda}$', 'Interpreter','latex')



% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])
%%  plot all mie outputs at 3 radii across all wavelengths

% plot real and imaginary refractive index
% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter
% plot the scattering efficiency

idx_radii = [1, 10, 50, 100];       % radii to plot

figure;
% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,1)
semilogy(ds.wavelength, ds.refrac_real)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Real part of Refractive Index', 'Interpreter','latex',...
    'FontSize', 20)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,2)
semilogy(ds.wavelength, ds.refrac_imag)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Imaginary part of Refractive Index', 'Interpreter','latex', ...
    'FontSize',20)



% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,3)
semilogy(ds.wavelength, ds.Qsca(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{sca}$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% in the second row, plot the extinction efficiency, the single scattering
% albedo and the asymmetry parameter
subplot(2,3,4)
semilogy(ds.wavelength, ds.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,5)
semilogy(ds.wavelength, ds.ssa(:,idx_radii))
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,6)
semilogy(ds.wavelength, ds.asymParam(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])


%% Plot just extinction coefficient, single scattering albedo and the asymmetry parameter


% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter

idx_radii = [1, 10, 50, 100];       % radii to plot


figure;
% plot the extinction efficiency
subplot(1,3,1)
semilogy(ds.wavelength, ds.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds.wavelength, ds.ssa(:,idx_radii))
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds.wavelength, ds.asymParam(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])



%% Plot extinction coefficient, single scattering albedo and the asymmetry parameter across all radii at a single wavelength


% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter

wavelen_2plot = 550;       % nm
wavelen_idx = ds.wavelength==wavelen_2plot;


figure;
% plot the extinction efficiency
subplot(1,3,1)
semilogy(ds.r_eff, ds.Qext(wavelen_idx,:), 'LineWidth',1, ...
    'Color', mySavedColors(1, 'fixed'))
grid on; grid minor
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex', 'Fontsize', 32)


% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds.r_eff, ds.ssa(wavelen_idx,:), 'LineWidth', 1, ...
    'Color', mySavedColors(1, 'fixed'))
grid on; grid minor
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex', 'Fontsize', 32)
title(['$\lambda = $',num2str(wavelen_2plot), ' $nm$'] ,'Interpreter','latex')



% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds.r_eff, ds.asymParam(wavelen_idx,:), 'LineWidth', 1, ...
    'Color', mySavedColors(1, 'fixed'))
grid on; grid minor
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex', 'Fontsize', 32)


% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])


%% Plot just extinction coefficient for a few effective radii


idx_radii = [1, 10, 50, 100];       % radii to plot


figure;
% plot the extinction efficiency
semilogy(ds.wavelength, ds.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)


% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])


%% Plot the extinction coefficient as a scaled image plot, showing all radii and all wavelengths



figure;
colormap hot
% plot the extinction efficiency
imagesc(ds.r_eff, ds.wavelength, ds.Qext)
colorbar
ylabel('Wavelength (nm)', 'Interpreter','latex')
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
title('$Q_{ext}$', 'Interpreter','latex')



% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])
