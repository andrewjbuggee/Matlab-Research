%% Lets compare the monodispersed distribution with the gamma droplet distribution

wl = 550;                       % nanometers
r = 1:100;                   % microns


% These calculations use the Wiscombe algorithm for computing Mie optical
% properties
mono = interp_mie_computed_tables([linspace(wl,wl,length(r))', r'], 'mono',false);

% These calculations use the Wiscombe algorithm for computing Mie optical
% properties
Gamma = interp_mie_computed_tables([linspace(wl,wl,length(r))', r'], 'gamma',false);

figure; subplot(1,3,1); semilogy(r,mono(:,5),r,Gamma(:,5),'LineWidth',2);
title('$Q_e$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); legend('monodispersed','gamma distribution')

subplot(1,3,2); plot(r,mono(:,6),r,Gamma(:,6),'LineWidth',2); 
title('$\tilde{\omega}$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$\tilde{\omega}$','Interpreter','latex');

subplot(1,3,3); plot(r,mono(:,7),r,Gamma(:,7),'LineWidth',2); 
title('$g$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$g$','Interpreter','latex')


%% Let's calculate Mie optical properties using the Wiscombe code for monodispersed
% and gamma dispersed distributions. Then let's compute an average over the
% droplet size distribution


wl = 550;                       % nanometers
r = 1:100;                   % microns

% These calculations use the Wiscombe algorithm for computing Mie optical
% properties
mono = interp_mie_computed_tables([linspace(wl,wl,length(r))', r'], 'mono',false);

% These calculations use the Wiscombe algorithm for computing Mie optical
% properties
Gamma = interp_mie_computed_tables([linspace(wl,wl,length(r))', r'], 'gamma',false);


% For some effective radius, we have defined a droplet size distribution
% ----- for a gamma distribution -----

r_modal = 1:100;                                              % microns
r = linspace(1, 100, 100);                  % microns - vector based on C.Emde (2016)
alpha = 7;
mu = alpha+3;                                            % to ensure I have the correct gamma distribution

ssa_avg = zeros(1,length(r_modal));
Qe_avg = zeros(1, length(r_modal));
Qe_avg2 = zeros(1, length(r_modal));

for ii = 1:length(r_modal)

    b = mu/r_modal(ii);                                   % exponent parameter
    N = mu^(mu+1)/(gamma(mu+1) * r_modal(ii)^(mu+1));  % normalization constant

    n_r = N*r.^mu .* exp(-b*r);                            % gamma droplet distribution



    % according to C.Emde (2016) page 1665, the average single scattering
    % albedo over a droplet distribution is:

    ssa_avg(ii) = trapz(r, pi*r.^2 .* mono(:,6)' .* mono(:,5)' .* n_r)./...
        trapz(r, pi*r.^2 .* mono(:,5)' .* n_r);

  
    Qe_avg(ii) = trapz(r, mono(:,8)' .* n_r)./trapz(r, mono(:,7)' .* n_r);

    Qe_avg2(ii) = trapz(r, mono(:,6)' .* n_r)./trapz(r, n_r);

end


figure; subplot(1,3,1); semilogy(r,mono(:,5),r,Gamma(:,5),'LineWidth',2);
hold on
semilogy(r,Qe_avg, r, Qe_avg2, 'LineWidth',2)
title('$Q_e$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); legend('monodispersed','gamma distribution')

subplot(1,3,2); plot(r,mono(:,6),r,Gamma(:,6),'LineWidth',2); 
hold on
semilogy(r,ssa_avg, 'LineWidth',2)
title('$\tilde{\omega}$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$\tilde{\omega}$','Interpreter','latex');

subplot(1,3,3); plot(r,mono(:,7),r,Gamma(:,7),'LineWidth',2); 
title('$g$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$g$','Interpreter','latex')


%% Lets run the Bohren and Huffman mie code for mono dispersed particles, gamma particles and log-normal particles


% What computer are you using? Grab the correct folder
if strcmp(whatComputer(),'anbu8374')==true
    folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
elseif strcmp(whatComputer(),'andrewbuggee')==true
    folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];
end

% Running bohren and Huffman for monodispersed
inputName = 'Mie_calcs_monodispersed_BH.INP';
outputName = 'OUTPUT_Mie_calcs_monodispersed_BH';
[drop_settings] = runMIE(folderName,inputName,outputName);
[BH_mono,~,~] = readMIE(folderName,outputName);



% Running bohren and Huffman for gamma distribution
inputName = 'Mie_calcs_gamma7_BH.INP';
outputName = 'OUTPUT_Mie_calcs_gamma7_BH';
[drop_settings] = runMIE(folderName,inputName,outputName);
[BH_gamma,~,~] = readMIE(folderName,outputName);


% Running bohren and Huffman for lognormal distribution
folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
inputName = 'Mie_calcs_logNormal_1_7_BH.INP';
outputName = 'OUTPUT_Mie_calcs_logNormal_1_7_BH';
[drop_settings] = runMIE(folderName,inputName,outputName);
[BH_logNorm,~,~] = readMIE(folderName,outputName);


%% lets make a plot comparing the Bohren-Huffman calculations for a single particle, a 
% gamma distribution and a log-normal  distribution


% What computer are you using? Grab the correct folder
if strcmp(whatComputer(),'anbu8374')==true
    folderName = '/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';
elseif strcmp(whatComputer(),'andrewbuggee')==true
    folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];
end

% Running bohren and Huffman for monodispersed
outputName = 'OUTPUT_Mie_calcs_monodispersed_BH';
[BH_mono,~,~] = readMIE(folderName,outputName);


% Running bohren and Huffman for gamma distribution
outputName = 'OUTPUT_Mie_calcs_gamma7_BH';
[BH_gamma,~,~] = readMIE(folderName,outputName);


% Running bohren and Huffman for lognormal distribution
outputName = 'OUTPUT_Mie_calcs_logNormal_1_7_BH';
[BH_logNorm,~,~] = readMIE(folderName,outputName);



r = 1:100;
wl = 550;   % nanometers

index = BH_mono.wavelength == wl;

figure; 
subplot(1,3,1); semilogy(r,BH_mono.Qext(index,:))
hold on; plot(r,BH_gamma.Qext(index,:),'linewidth',2)
plot(r,BH_logNorm.Qext(index,:), 'linewidth',2);
title('$Q_e$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); legend('monodispersed','gamma','log normal')

subplot(1,3,2); plot(r,BH_mono.ssa(index,:))
hold on; plot(r,BH_gamma.ssa(index,:),'linewidth',2)
plot(r,BH_logNorm.ssa(index,:),'linewidth',2); 
title('$\tilde{\omega}$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$\tilde{\omega}$','Interpreter','latex');

subplot(1,3,3); plot(r,BH_mono.asymParam(index,:))
hold on; plot(r,BH_gamma.asymParam(index,:),'linewidth',2)
plot(r,BH_logNorm.asymParam(index,:),'linewidth',2);
title('$g$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$g$','Interpreter','latex')





%% Attempting to solve for the single scattering albedo of a droplet distribution

% define the wavelength of interest
wl = 550;   % nanometers

index = BH_mono.wavelength == wl;

% For some effective radius, we have defined a droplet size distribution
% ----- for a gamma distribution -----

% r needs to span 5 orders of magnitude to ensure the normaliztion N
% will cause the integral to equal 1 for all modal radii used
r = logspace(-2,3,1000);
r_modal = 1:100;                                              % microns
alpha = 7;                                              % alpha = 7 is libRadTran's cited number to use for liquid water droplets
mu = alpha+3;                                            % to ensure I have the correct gamma distribution

ssa_avg = zeros(1,length(r_modal));
ssa_avg2 = zeros(1, length(r_modal));
Qe_avg = zeros(1, length(r_modal));
Qe_avg2 = zeros(1, length(r_modal));

r_eff_calc = zeros(1, length(r_modal));
r_modal_calc = zeros(1, length(r_modal));
Nc_calc = zeros(1, length(r_modal));


for ii = 1:length(r_modal)

    b = mu/r_modal(ii);                                   % exponent parameter
    N = mu^(mu+1)/(gamma(mu+1) * r_modal(ii)^(mu+1));     % normalization constant
    


    % This distribution integrates to 1?
    n_r = N*r.^mu .* exp(-b*r);                            % gamma droplet distribution

    % do I get the effective radius that I used to define the distribution?
    r_eff_calc(ii) = trapz(r, r.^3 .* n_r)./trapz(r, r.^2 .* n_r);

    % what's the modal radius for each distribution?
    [~,modal_index] = max(n_r);
    r_modal_calc(ii) = r(modal_index);

    % does my distribution integrate to 1?
    Nc_calc(ii) = trapz(r, n_r);

    % according to C.Emde (2016) page 1665, the average single scattering
    % albedo over a droplet distribution is:

    ssa_avg(ii) = trapz(r, (pi*r.^2) .* (BH_mono.ssa(index,:)) .* (BH_mono.Qext(index,:) .* (n_r)))./...
        trapz(r, (pi*r.^2) .* (BH_mono.Qext(index,:) .* (n_r)));

    ssa_avg2(ii) = trapz(r, BH_mono.ssa(index,:) .* n_r)./trapz(r, n_r);

    Qe_avg(ii) = trapz(r, (BH_mono.Qsca(index,:)) .* n_r)./trapz(r, BH_mono.ssa(index,:) .* n_r);

    Qe_avg2(ii) = trapz(r, BH_mono.Qext(index,:) .* n_r)./trapz(r, n_r);


end




%% Now make a plot comparing my computations of Q_e and omega for a gamma droplet distribution with libRadTran

% So let's only plot libRadTran computations for a single particle and for
% a gamma distribution

% Grab th seafoam green from mySavedColors
C_gold = mySavedColors(7, 'fixed');

nonMonoWidth = 3;

index = BH_mono.wavelength == wl;

figure; 
subplot(1,3,1); semilogy(r,BH_mono.Qext(index,:))
hold on; plot(r,BH_gamma.Qext(index,:),'linewidth',nonMonoWidth)
plot(r,Qe_avg, '--', 'linewidth',nonMonoWidth,'Color',C_gold);
title('$Q_e$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); 
legend('monodispersed','gamma - libradTran','gamma - my calculation')

subplot(1,3,2); plot(r,BH_mono.ssa(index,:))
hold on; plot(r,BH_gamma.ssa(index,:),'linewidth',nonMonoWidth)
plot(r,ssa_avg, '--', 'linewidth',nonMonoWidth,'Color',C_gold);
plot(r, ssa_avg2, 'k--', 'LineWidth',nonMonoWidth)
title('$\tilde{\omega}$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$\tilde{\omega}$','Interpreter','latex');
legend('monodispersed','gamma - libradTran','gamma - Emde avg', 'gamma - simple avg')


subplot(1,3,3); plot(r,BH_mono.asymParam(index,:))
hold on; plot(r,BH_gamma.asymParam(index,:),'linewidth',nonMonoWidth)
title('$g$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$g$','Interpreter','latex')

%% Plot only libRadTran's single particle mie calculation for the extinction efficiency, along with my estimate 
% for the extinction efficiency of a gamma droplet distribution


%% Now make a plot comparing my computations of Q_e and omega for a gamma droplet distribution with libRadTran

% So let's only plot libRadTran computations for a single particle and for
% a gamma distribution

% Grab th seafoam green from mySavedColors
C_gold = mySavedColors(7, 'fixed');
C_pink = mySavedColors(1, 'fixed');

nonMonoWidth = 3;

index = BH_mono.wavelength == wl;

figure; 
plot(r,BH_mono.Qext(index,:))
hold on; 
plot(r,Qe_avg, '--', 'linewidth',nonMonoWidth,'Color',C_gold);
title('$Q_e$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); 
legend('monodispersed','gamma - my calculation')



figure; 
plot(r,BH_mono.Qext(index,:))
hold on; 
plot(r,Qe_avg, '--', 'linewidth',nonMonoWidth,'Color',C_gold);
plot(r,Qe_avg2, '--', 'linewidth',nonMonoWidth,'Color',C_pink);
title('$Q_e$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); 
legend('monodispersed','gamma - calc1', 'gamma - calc2')


%% Let's compare Qe and Qscat from libRadTran

index_abs = BH_mono.wavelength==1900;

figure; 
subplot(1,2,1); semilogy(r,BH_mono.Qext(index_abs,:))
hold on; plot(r,BH_gamma.Qext(index_abs,:),'linewidth',nonMonoWidth)
plot(r,BH_logNorm.Qext(index_abs,:), 'linewidth',nonMonoWidth);
title('$Q_e$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_e$','Interpreter','latex'); 
legend('monodispersed','gamma - libradTran','log-noraml - libRadTran')

subplot(1,2,2); semilogy(r,BH_mono.Qsca(index_abs,:))
hold on; plot(r,BH_gamma.Qsca(index_abs,:),'linewidth',nonMonoWidth)
plot(r,BH_logNorm.Qsca(index_abs,:), 'linewidth',nonMonoWidth);
title('$Q_{s}$ using BH code','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$)','Interpreter','latex'); ylabel('$Q_s$','Interpreter','latex');


%% Plot my calculated effective radius and total number concentration

figure; 
subplot(1,2,1); plot(r,r_eff_calc)
hold on; plot(r, r_modal_calc)
plot(r,r,'k--')
title('Calculated $r_{eff}$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$) used to define $n(r)$','Interpreter','latex'); 
ylabel('Calculated $r_e$','Interpreter','latex'); 
legend('$r_e$', 'modal $r$', '$1:1$','Interpreter','latex', 'location','best'); 


subplot(1,2,2); plot(r,Nc_calc)
title('Calculated $N_{c}$','Interpreter','latex'); grid on; grid minor
xlabel('$r_{eff}$ ($\mu m$) used to define $n(r)$','Interpreter','latex'); 
ylabel('$N_c$','Interpreter','latex'); 
