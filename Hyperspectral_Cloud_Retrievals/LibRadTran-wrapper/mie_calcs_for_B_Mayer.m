%% Plots for my Email to Berhard Mayer


% Andrew John Bugee

%%  When using LibRadTran's Mie solver, if one assumes a droplet distribution, the calcualtions grow extremly large

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
ylabel('$\curlypi_0$', 'Interpreter','latex')
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
