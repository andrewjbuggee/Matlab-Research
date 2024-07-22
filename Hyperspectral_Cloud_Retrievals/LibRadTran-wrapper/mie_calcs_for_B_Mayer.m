%% Plots for my Email to Berhard Mayer


% Andrew John Bugee

%%  When using LibRadTran's Mie solver, if one assumes a droplet distribution, the calcualtions of extinction efficiency
% tend towards infinity as the effective radius tends towards 0

clear variables

%% ----- Run Monodispersed Mie File using MIEV0 solver -----


folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_B-Mayer_mono_dispersed_MIEV0.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_miev0,headers,num_radii] = readMIE(folderName,inputName(1:end-4));



%%  plot all mie outputs all radii across 5 wavelengths
% This is for a file that ran with only 5 wavelengths

% plot real and imaginary refractive index
% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter
% plot the scattering efficiency


figure;
% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,1)
semilogy(ds_miev0.wavelength, ds_miev0.refrac_real)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Real part of Refractive Index', 'Interpreter','latex',...
    'FontSize', 20)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,2)
semilogy(ds_miev0.wavelength, ds_miev0.refrac_imag)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Imaginary part of Refractive Index', 'Interpreter','latex', ...
    'FontSize',20)



% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,3)
semilogy(ds_miev0.r_eff, ds_miev0.Qsca, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{sca}$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% in the second row, plot the extinction efficiency, the single scattering
% albedo and the asymmetry parameter
subplot(2,3,4)
semilogy(ds_miev0.r_eff, ds_miev0.Qext, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,5)
semilogy(ds_miev0.r_eff, ds_miev0.ssa)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)


% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,6)
semilogy(ds_miev0.r_eff, ds_miev0.asymParam, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])



%%  plot 3 mie outputs all radii across 5 wavelengths for monodispersed droplets
% This is for a file that ran with only 5 wavelengths

% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter


figure;

% plot the extinction efficiency
subplot(1,3,1)
semilogy(ds_miev0.r_eff, ds_miev0.Qext, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds_miev0.r_eff, ds_miev0.ssa, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)
% define plot title
title('Monodispersed Calcuations using MIEV0', 'Interpreter','latex')


% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds_miev0.r_eff, ds_miev0.asymParam, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_miev0.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])





%% ----- Run Monodispersed Mie File using Bohren-Huffman solver -----


folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_B-Mayer_mono_dispersed_BH.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_bh,headers,num_radii] = readMIE(folderName,inputName(1:end-4));



%%  plot 3 mie outputs all radii across 5 wavelengths for monodispersed droplets
% This is for a file that ran with only 5 wavelengths

% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter


figure;

% plot the extinction efficiency
subplot(1,3,1)
semilogy(ds_bh.r_eff, ds_bh.Qext, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_bh.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds_bh.r_eff, ds_bh.ssa, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_bh.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)
% define plot title
title('Monodispersed Calcuations using MIEV0', 'Interpreter','latex')


% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds_bh.r_eff, ds_bh.asymParam, 'LineWidth',2)
grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(append('$\lambda = $',string(ds_bh.wavelength), ' $nm$'), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])







%% ----- Run Gamma and lognormal distribution Mie Files using MIEV0 solver -----


folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_B-Mayer_gamma_MIEV0.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_miev0_gamma,headers,num_radii] = readMIE(folderName,inputName(1:end-4));


inputName = 'Mie_calcs_for_B-Mayer_lognormal_MIEV0.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_miev0_lognormal,headers,num_radii] = readMIE(folderName,inputName(1:end-4));







%% Compare MIEV0 mie outputs for all radii across 5 wavelengths for monodispersed and gamma dispersed droplets
% This is for a file that ran with only 5 wavelengths
% --- Plot each seperately ---

% plot the extinction efficiency
% plot the scattering efficiency
% plot the absorption efficiency
% plot the single scattering albedo
% plot the asymmetry parameter




% --- first plot gamma dispersed droplets, then monodispersed ---

% plot the extinction efficiency
legend_str = cell(3*length(ds_miev0.wavelength), 1);

figure;

for nn = 1:length(ds_miev0.wavelength)

    semilogy(ds_miev0.r_eff, ds_miev0.Qext(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_miev0_gamma.r_eff, ds_miev0_gamma.Qext(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    semilogy(ds_miev0_lognormal.r_eff, ds_miev0_lognormal.Qext(nn,:), 'LineStyle', '--', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{3*nn -2} = ['Mono; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn -1} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn} = ['Log Normal; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Extinction Efficiency using MIEV0 - $\alpha=7$', 'Interpreter','latex')
% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the single scattering albedo
legend_str = cell(3*length(ds_miev0.wavelength), 1);

figure;

for nn = 1:length(ds_miev0.wavelength)

    semilogy(ds_miev0.r_eff, ds_miev0.ssa(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_miev0_gamma.r_eff, ds_miev0_gamma.ssa(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    semilogy(ds_miev0_lognormal.r_eff, ds_miev0_lognormal.ssa(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{3*nn -2} = ['Mono; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn -1} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn} = ['Log Normal; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$\varpi$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Single Scattering Albedo using MIEV0 - $\alpha=7$', 'Interpreter','latex')
% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])



% plot the asymmetry parameter
legend_str = cell(3*length(ds_miev0.wavelength), 1);

figure;

for nn = 1:length(ds_miev0.wavelength)

    semilogy(ds_miev0.r_eff, ds_miev0.asymParam(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_miev0_gamma.r_eff, ds_miev0_gamma.asymParam(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

     semilogy(ds_miev0_lognormal.r_eff, ds_miev0_lognormal.asymParam(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{3*nn -2} = ['Mono; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn -1} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn} = ['Log Normal; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Asymmetry Parameter using MIEV0 - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the scattering efficiency
legend_str = cell(3*length(ds_miev0.wavelength), 1);

figure;

for nn = 1:length(ds_miev0.wavelength)

    semilogy(ds_miev0.r_eff, ds_miev0.Qsca(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_miev0_gamma.r_eff, ds_miev0_gamma.Qsca(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    semilogy(ds_miev0_lognormal.r_eff, ds_miev0_lognormal.Qsca(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));


    legend_str{3*nn -2} = ['Mono; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn -1} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn} = ['Log Normal; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{sca}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Scattering Efficiency using MIEV0 - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the absorption efficiency
legend_str = cell(3*length(ds_miev0.wavelength), 1);

figure;

for nn = 1:length(ds_miev0.wavelength)

    semilogy(ds_miev0.r_eff, ds_miev0.Qext(nn,:) - ds_miev0.Qsca(nn,:), 'LineStyle', ':',...
        'LineWidth',2,'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_miev0_gamma.r_eff, ds_miev0_gamma.Qext(nn,:) - ds_miev0_gamma.Qsca(nn,:),...
        'LineStyle', '-', 'LineWidth',2,'Color', mySavedColors(nn, 'fixed'));

    semilogy(ds_miev0_lognormal.r_eff, ds_miev0_lognormal.Qext(nn,:) - ds_miev0_gamma.Qsca(nn,:),...
        'LineStyle', '-', 'LineWidth',2,'Color', mySavedColors(nn, 'fixed'));

    legend_str{3*nn -2} = ['Mono; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn -1} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];
    legend_str{3*nn} = ['Gamma; $\lambda = $',num2str(ds_miev0_gamma.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{abs}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Absorption Efficiency using MIEV0 - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])






%% ----- Run Mono and Gamma dispersed Mie Files using Bohren-Huffman solver -----


folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_B-Mayer_mono_dispersed_BH.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_bh,headers,num_radii] = readMIE(folderName,inputName(1:end-4));


inputName = 'Mie_calcs_for_B-Mayer_gamma_BH.INP';

tic
runMIE(folderName,inputName,inputName(1:end-4));
toc

% read .OUT file
[ds_bh_gamma,headers,num_radii] = readMIE(folderName,inputName(1:end-4));





%% Compare Bohren-Huffman mie outputs for all radii across 5 wavelengths for monodispersed and gamma dispersed droplets
% This is for a file that ran with only 5 wavelengths
% --- Plot each seperately ---

% plot the extinction efficiency
% plot the scattering efficiency
% plot the absorption efficiency
% plot the single scattering albedo
% plot the asymmetry parameter




% --- first plot gamma dispersed droplets, then monodispersed ---

% plot the extinction efficiency
legend_str = cell(2*length(ds_bh.wavelength), 1);

figure;

for nn = 1:length(ds_bh.wavelength)

    semilogy(ds_bh.r_eff, ds_bh.Qext(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_bh_gamma.r_eff, ds_bh_gamma.Qext(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{2*nn -1} = ['Mono; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];
    legend_str{2*nn} = ['Gamma; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Extinction Efficiency using Bohren-Huffman - $\alpha=7$', 'Interpreter','latex')
% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the single scattering albedo
legend_str = cell(2*length(ds_bh.wavelength), 1);

figure;

for nn = 1:length(ds_bh.wavelength)

    semilogy(ds_bh.r_eff, ds_bh.ssa(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_bh_gamma.r_eff, ds_bh_gamma.ssa(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{2*nn -1} = ['Mono; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];
    legend_str{2*nn} = ['Gamma; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$\varpi$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Single Scattering Albedo using Bohren-Huffman - $\alpha=7$', 'Interpreter','latex')
% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])



% plot the asymmetry parameter
legend_str = cell(2*length(ds_bh.wavelength), 1);

figure;

for nn = 1:length(ds_bh.wavelength)

    semilogy(ds_bh.r_eff, ds_bh.asymParam(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_bh_gamma.r_eff, ds_bh_gamma.asymParam(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{2*nn -1} = ['Mono; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];
    legend_str{2*nn} = ['Gamma; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Asymmetry Parameter using Bohren-Huffman - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the scattering efficiency
legend_str = cell(2*length(ds_bh.wavelength), 1);

figure;

for nn = 1:length(ds_bh.wavelength)

    semilogy(ds_bh.r_eff, ds_bh.Qsca(nn,:), 'LineStyle', ':', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_bh_gamma.r_eff, ds_bh_gamma.Qsca(nn,:), 'LineStyle', '-', 'LineWidth',2,...
        'Color', mySavedColors(nn, 'fixed'));

    legend_str{2*nn -1} = ['Mono; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];
    legend_str{2*nn} = ['Gamma; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{sca}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Scattering Efficiency using Bohren-Huffman - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




% plot the absorption efficiency
legend_str = cell(2*length(ds_bh.wavelength), 1);

figure;

for nn = 1:length(ds_bh.wavelength)

    semilogy(ds_bh.r_eff, ds_bh.Qext(nn,:) - ds_bh.Qsca(nn,:), 'LineStyle', ':',...
        'LineWidth',2,'Color', mySavedColors(nn, 'fixed'));
    hold on
    semilogy(ds_bh_gamma.r_eff, ds_bh_gamma.Qext(nn,:) - ds_bh_gamma.Qsca(nn,:),...
        'LineStyle', '-', 'LineWidth',2,'Color', mySavedColors(nn, 'fixed'));

    legend_str{2*nn -1} = ['Mono; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];
    legend_str{2*nn} = ['Gamma; $\lambda = $',num2str(ds_bh.wavelength(nn)), ' $nm$'];

end

grid on; grid minor
xlabel('Radius ($\mu m$)', 'Interpreter','latex')
ylabel('$Q_{abs}$', 'Interpreter','latex')
% legend([append('Gamma; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$'); ...
%    append('Mono; $\lambda = $',string(ds_miev0_gamma.wavelength), ' $nm$')] , 'Location',...
%    'best','Interpreter','latex', 'FontSize',15)
legend(legend_str, 'Location',...
   'best','Interpreter','latex', 'FontSize',15)
% define plot title
title('Absorption Efficiency using Bohren-Huffman - $\alpha=7$', 'Interpreter','latex')

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])




%%





%%  Plot the size parameter across all radii and wavelength values

figure; 
imagesc(ds_miev0.r_eff, ds_miev0.wavelength, 2*pi*ds_miev0.r_eff./(ds_miev0.wavelength./1e3));
colormap hot
% plot the extinction efficiency
colorbar
ylabel('Wavelength (nm)', 'Interpreter','latex')
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
title('Size Parameter - $x = \frac{2 \pi r}{\lambda}$', 'Interpreter','latex')



% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])


%%  plot all mie outputs at 5 radii across all wavelengths

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
semilogy(ds_miev0.wavelength, ds_miev0.refrac_real)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Real part of Refractive Index', 'Interpreter','latex',...
    'FontSize', 20)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,2)
semilogy(ds_miev0.wavelength, ds_miev0.refrac_imag)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
title('Imaginary part of Refractive Index', 'Interpreter','latex', ...
    'FontSize',20)



% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,3)
semilogy(ds_miev0.wavelength, ds_miev0.Qsca(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{sca}$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% in the second row, plot the extinction efficiency, the single scattering
% albedo and the asymmetry parameter
subplot(2,3,4)
semilogy(ds_miev0.wavelength, ds_miev0.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,5)
semilogy(ds_miev0.wavelength, ds_miev0.ssa(:,idx_radii))
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% in the first row, plot the refractive indices, and the scattering
% efficiency. 
subplot(2,3,6)
semilogy(ds_miev0.wavelength, ds_miev0.asymParam(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
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
semilogy(ds_miev0.wavelength, ds_miev0.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds_miev0.wavelength, ds_miev0.ssa(:,idx_radii))
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)



% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds_miev0.wavelength, ds_miev0.asymParam(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$g$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)

% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])



%% Plot extinction coefficient, single scattering albedo and the asymmetry parameter across all radii at a single wavelength


% plot the extinction efficiency
% plot the single scattering albedo
% plot the asymmetry parameter

wavelen_2plot = 550;       % nm
wavelen_idx = ds_miev0.wavelength==wavelen_2plot;


figure;
% plot the extinction efficiency
subplot(1,3,1)
semilogy(ds_miev0.r_eff, ds_miev0.Qext(wavelen_idx,:), 'LineWidth',1, ...
    'Color', mySavedColors(1, 'fixed'))
grid on; grid minor
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex', 'Fontsize', 32)


% plot the single scattering albedo
subplot(1,3,2)
semilogy(ds_miev0.r_eff, ds_miev0.ssa(wavelen_idx,:), 'LineWidth', 1, ...
    'Color', mySavedColors(1, 'fixed'))
grid on; grid minor
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
ylabel('$\varpi_0$', 'Interpreter','latex', 'Fontsize', 32)
title(['$\lambda = $',num2str(wavelen_2plot), ' $nm$'] ,'Interpreter','latex')



% plot the asymmetry parameter
subplot(1,3,3)
semilogy(ds_miev0.r_eff, ds_miev0.asymParam(wavelen_idx,:), 'LineWidth', 1, ...
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
semilogy(ds_miev0.wavelength, ds_miev0.Qext(:,idx_radii), 'LineWidth',2)
grid on; grid minor
xlabel('Wavelength (nm)', 'Interpreter','latex')
ylabel('$Q_{ext}$', 'Interpreter','latex')
legend(string(ds_miev0.r_eff(idx_radii)), 'Location','best','Interpreter','latex',...
    'FontSize',15)


% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])


%% Plot the extinction coefficient as a scaled image plot, showing all radii and all wavelengths



figure;
colormap hot
% plot the extinction efficiency
imagesc(ds_miev0.r_eff, ds_miev0.wavelength, ds_miev0.Qext)
colorbar
ylabel('Wavelength (nm)', 'Interpreter','latex')
xlabel('Effective Radius $(\mu m)$', 'Interpreter','latex')
title('$Q_{ext}$', 'Interpreter','latex')



% Set general figure settings
set(gcf, 'Position', [0 0 1250, 600])
