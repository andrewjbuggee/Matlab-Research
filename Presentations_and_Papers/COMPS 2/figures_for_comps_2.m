% FIGURES MADE FOR COMPS 2


%% Plot droplet distribution at the base of some vert profile

clear variables

% Load some VOCALS-REx Data

% Oct-18-2008 Data
filename = 'RF02.20081018.130300_213000.PNI.nc';

% Load data
vocalsRex = readVocalsRex(filename);


% Plot one distribuiton at the bottom of the cloud
vocals_time1 = 5916-1;
index_time1 = vocalsRex.time==vocals_time1;


% Compute the effective radius and plot it as a solid vertical line
re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
      sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');

% Compute total droplet number concentration
Nc1 = sum(vocalsRex.Nc(:,index_time1));


% plot another distribution at the top of the cloud
vocals_time2 = 5964-1;
index_time2 = vocalsRex.time==vocals_time2;



% Compute the effective radius and plot it as a solid vertical line
re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
      sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');

% Compute total droplet number concentration
Nc2 = sum(vocalsRex.Nc(:,index_time2));

% Plot both of them on the same figure!


f1 = figure;
h1 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1)); 
h1.FaceColor = mySavedColors(11, 'fixed');
h1.FaceAlpha = 0.7;
h1.EdgeAlpha = 1;
hold on
xline(re1,'k--', 'LineWidth',4)

h2 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2)); 
h2.FaceColor = mySavedColors(14, 'fixed');
h2.FaceAlpha = 0.4;
h2.EdgeAlpha = 1;
xline(re2,'k--', 'LineWidth',4)

xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
ylabel('$N$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize',32);
grid on; grid minor; hold on;
xlim([0,25])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
set(gcf, 'Position',[0 0 750, 450])

% Label first xline
annotation(f1,'textbox',[0.166458333333333 0.92 0.0911875000000001 0.0711111111111111],...
    'String','$r_e = 4.7 \mu m$',...
    'Interpreter','latex',...
    'Color',mySavedColors(11,'fixed'),...
    'FontSize',24,...
    'FontWeight','bold',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label second x line
annotation(f1,'textbox',[0.340072916666669 0.92 0.0911875000000001 0.0711111111111111],...
    'String','$r_e = 7.2 \mu m$',...
    'Interpreter','latex',...
    'Color',mySavedColors(4,'fixed'),...
    'FontSize',24,...
    'FontWeight','bold',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label the first number concentration
annotation(f1,'textbox',[0.642739583333335 0.737777777777777 0.225260416666665 0.0711111111111111],...
    'String',['$N_c = $',num2str(Nc1,3),' $cm^{-3}$'],...
    'Interpreter','latex',...
    'Color',mySavedColors(11,'fixed'),...
    'FontSize',24,...
    'FontWeight','bold',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label the second number concentration
annotation(f1,'textbox',[0.642739583333335 0.817777777777778 0.230593749999998 0.0711111111111111],...
    'String',['$N_c = $',num2str(Nc2,3),' $cm^{-3}$'],...
    'Interpreter','latex',...
    'Color',mySavedColors(4,'fixed'),...
    'FontSize',24,...
    'FontWeight','bold',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(gca,'FontSize',20)



%% Create plots of simple Monte Carlo demonstrations

clear variables


% define the original distribution

P_tau = @(tau) exp(-tau);

% Define the uniform distribution

N_samples1 = 100000;
N_samples2 = 1000;

% reset the random number generator
rng('default');


% define the relationship between x and tau

tau = @(x) -log(1 - x);


% We need to define the bins for the distribution we wish to compute. But
% the range is semi-infinite. So how do we treat this? We will allow the
% last bin to go to infinity. Usually, for N bins we have to define N+1
% edge points. In this case, we will only define N edge points and allow
% the last edge point to go to infinity

tau_boundaries = [0, inf];

% our bins are uniformly distributed
bin_width = 0.05; 

% define the bin edges
bin_edges = [tau_boundaries(1):bin_width:10, tau_boundaries(2)];

% define the number of bins

N_bins = length(bin_edges)-1;

% define our new distribution
N_counts1 = zeros(1,N_bins);

% assign each x to one of our bins
for nn=1:N_samples1

    % draw a value from our uniform distribution
    x_draw = rand(1,1);
    tau_sample = tau(x_draw);

    for bb = 1:(N_bins+1)

        if tau_sample>=bin_edges(bb) && tau_sample<bin_edges(bb+1)

            N_counts1(bb) = N_counts1(bb) + 1;

        end
    end
end

PDF1 = (N_counts1./N_samples1)./bin_width;


% define our new distribution
N_counts2 = zeros(1,N_bins);
% assign each x to one of our bins
for nn=1:N_samples2

    % draw a value from our uniform distribution
    x_draw = rand(1,1);
    tau_sample = tau(x_draw);

    for bb = 1:(N_bins+1)

        if tau_sample>=bin_edges(bb) && tau_sample<bin_edges(bb+1)

            N_counts2(bb) = N_counts2(bb) + 1;

        end
    end
end

PDF2 = (N_counts2./N_samples2)./bin_width;


% Let's divide the number of counts in each bin by the total number of
% counts and the bin width to obtain the PDF

% Let's define the vector x to plot the PDF the center point of each bin
tau_PDF = bin_edges(1)+(bin_width/2):bin_width:bin_edges(end-1)+(bin_width/2);




% Let's plot this PDF on top of our original PDF to see how we did


figure; plot(tau_PDF, P_tau(tau_PDF),"Color",mySavedColors(13,'fixed'),...
    'LineWidth',6)
hold on;
plot(tau_PDF, PDF2,"Color",mySavedColors(15,'fixed'),'LineWidth',3)
plot(tau_PDF, PDF1,"Color",'k','LineWidth',3)
grid on; grid minor
xlabel('$\tau$','Interpreter','latex', 'FontSize',32);
ylabel('$P(\tau)$','Interpreter','latex', 'FontSize',32);
legend('$P(\tau)$',['\# Samples: ',num2str(N_samples2)]',['\# Samples: ',num2str(N_samples1)]',...
    'interpreter','latex','Location','best','FontSize',24)



set(gcf, 'Position',[0 0 750 350])
set(gca,'FontSize',20)


%% Plot of the bulk absorption coefficient of liquid water

clear variables; 

% load mie calc file with 1nm resolution across radii from 1 to 100 microns

folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_AVIRIS_monodispersed_BH.INP';

outputName = 'OUTPUT_Mie_calcs_for_AVIRIS_monodispersed_BH';

% read mie_calculation data
[ds,headers,num_radii] = readMIE(folderName,outputName);

% What wavelengths do you want to plot? Data ranges from 100-3000 nm
wavelength_2Plot = ds.wavelength>=350 & ds.wavelength<=2300;


% Make plot of the bulk absorption coefficient

figure;
semilogy(ds.wavelength(wavelength_2Plot), (4*pi*ds.refrac_imag(wavelength_2Plot))./(ds.wavelength(wavelength_2Plot)*1e-7),...
    "Color",mySavedColors(11,'fixed'),'LineWidth',4)
grid on; grid minor
xlabel('Wavelength ($nm$)','Interpreter','latex', 'FontSize',32);
ylabel('$\kappa$ ($cm^{-1}$)','Interpreter','latex', 'FontSize',32);
set(gcf, 'Position',[0 0 750 475])
set(gca,'FontSize',20)



%% Plot 1-ssa for 10 different droplet sizes

clear variables; 

% load mie calc file with 1nm resolution across radii from 1 to 100 microns

folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_AVIRIS_monodispersed_BH.INP';

outputName = 'OUTPUT_Mie_calcs_for_AVIRIS_monodispersed_BH';

% read mie_calculation data
[ds,headers,num_radii] = readMIE(folderName,outputName);


% Make plot of the bulk absorption coefficient
radii_2_plot = [1,10:10:100];
wavelength_2Plot = ds.wavelength<2300;


figure;

for ii = 1:length(radii_2_plot)

    radius_index = ds.r_eff==radii_2_plot(ii);

    plot(ds.wavelength(wavelength_2Plot), (1 - ds.ssa(wavelength_2Plot,radius_index)),...
     "Color",mySavedColors(ii,'fixed'),'LineWidth',2.5)

    hold on

    legend_str{ii} = ['$r_e = $',num2str(radii_2_plot(ii)), ' $\mu m$'];

end


grid on; grid minor
xlabel('Wavelength ($nm$)','Interpreter','latex', 'FontSize',32);
ylabel('$1 - \varpi_0$','Interpreter','latex', 'FontSize',32);
legend(legend_str,'interpreter','latex','Location','best','FontSize',20,...
    'Position',[0.20440577301917 0.301406822806967 0.174991692258035 0.613559325985091])
set(gcf, 'Position',[0 0 850 550])
set(gca,'FontSize',20)




%% Plot Q_scat and Q_abs for a single droplet size of 20 microns

clear variables; 

% load mie calc file with 1nm resolution across radii from 1 to 100 microns

folderName = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

inputName = 'Mie_calcs_for_AVIRIS_monodispersed_BH.INP';

outputName = 'OUTPUT_Mie_calcs_for_AVIRIS_monodispersed_BH';

% read mie_calculation data
[ds,headers,num_radii] = readMIE(folderName,outputName);


% Make plot of the bulk absorption coefficient
radii_2_plot = 20;
radius_index = ds.r_eff==radii_2_plot;
wavelength_2Plot = ds.wavelength>200 & ds.wavelength<2300;

figure; 
% first plot the scattering efficiency
ax1 = subplot(2,1,1);
semilogy(ds.wavelength(wavelength_2Plot)./1e3, ds.Qsca(wavelength_2Plot,radius_index),...
     "Color",mySavedColors(9,'fixed'),'LineWidth',3.5)
grid on; grid minor
ylabel('$Q_s$','Interpreter','latex', 'FontSize',32);
title(['Liquid water droplet: $r = $',num2str(radii_2_plot), ' $\mu m$'], ...
    'Interpreter','latex')
set(ax1,'Position',[0.13 0.555946295772496 0.775 0.316417346139791])


% Second, plot the absorption efficiency
ax2 = subplot(2,1,2);
semilogy(ds.wavelength(wavelength_2Plot)./1e3, (ds.Qext(wavelength_2Plot,radius_index) - ds.Qsca(wavelength_2Plot,radius_index)),...
        "Color",mySavedColors(9,'fixed'),'LineWidth',3.5)
grid on; grid minor
ylabel('$Q_a$','Interpreter','latex', 'FontSize',32);
xlabel('Wavelength ($\mu m$)','Interpreter','latex', 'FontSize',32);

set(gcf, 'Position',[0 0 850 550])
set(ax1,'FontSize',18)
set(ax2,'FontSize',18)
set(ax2,'YTick', [10^(-6), 10^(-4), 10^(-2), 10^(0)])
set(ax2,'YTickLabel', {'10^{-6}','10^{-4}','10^{-2}', '10^{0}'})



%% Compute several approximations to the fractional absorption of an infintely thick cloud
% Twomey and Bohren 1980

% Lets start by picking a representative solar zenith angle
theta0 = 0;
mu0 = cosd(theta0);         % cosine of the solar zenith angle


% Pick a wavelength range
wl_nm = 350:2300;                  % nm
wl_microns = wl_nm./1e3;            % microns

% Lets vary the droplet radius across a reasonable range that could be
% measured in a warm liquid cloud
r = 5;                     % microns

% compute the single scattering albedo
mie_prop = zeros(length(wl_nm), 8, length(r));
single_scattering_albedo = zeros(length(wl_nm), length(r));


g = zeros(length(wl_nm), length(r));

% We also need to compute the H function for all single scattering values
H_func = zeros(length(wl_nm), length(r));

for rr = 1:length(r)
    xq = [wl_nm',linspace(r(rr),r(rr),length(wl_nm))'];
    mie_prop(:,:,rr) = interp_mie_computed_tables(xq,'mono',false);
    
    % Grab the single scattering albedo
    single_scattering_albedo(:,rr) = reshape(mie_prop(:,6,rr),length(wl_nm),1);
    g(:,rr) = reshape(mie_prop(:,7,rr),length(wl_nm),1);

    % Smooth out single scattering albedo. We cannot measure the high frequency
    % content
    single_scattering_albedo(:,rr) = smoothdata(single_scattering_albedo(:,rr), 'gaussian',50);
    g(:,rr) = smoothdata(g(:,rr), 'gaussian', 50);

    % Compute the H function value
    H_func(:,rr) = H(mu0, single_scattering_albedo(:,rr), 0);
end


% Calculate fractional absorption using eq. 3, which doesn't use many
% simplifications
A = H_func.* sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo));

% Calculate the simpler absorption
% read in and interpolate values of the complex index of refraction
refrac_ind = interp_refractive_index_H2O(wl_microns);

% Compute the absorption coefficient
k = 4*pi*refrac_ind(:,2)./wl_microns';               % microns^(-1) - absorption coefficient of water

A_simp = H_func.* sqrt((0.85*k*r)./(1 - g.*(1 - 0.85*k*r)));

A_simp2 = 2.2*H_func.*sqrt(k*r);



% Make a plot!
figure; 
plot(wl_nm, A,'LineWidth',4, 'Color', mySavedColors(8,'fixed'))
hold on;
plot(wl_nm, A_simp,'LineWidth',4, 'Color', mySavedColors(9,'fixed'))
plot(wl_nm, A_simp2,'LineWidth',4, 'Color', mySavedColors(10,'fixed'))

grid on; grid minor
xlabel('Wavelength ($nm$)','Interpreter','latex')
ylabel('Fractional Absorption', 'Interpreter','latex')

legend({'$H(\mu_0) \sqrt{\frac{1 - \omega_0}{1 - g\omega_0}}$',...
    '$H(\mu_0) \sqrt{\frac{0.85 r\kappa}{1 - g(1 - 0.85 r \kappa)}}$',...
    '$2.2 H(\mu_0) \sqrt{r \kappa}$'},...
    'Interpreter','latex', 'FontSize',25,...
    'Position',[0.157647058823529 0.601818181818182 0.316811272116267 0.303609756007339]);


set(gcf, 'Position',[0 0 850 550])
set(gca,'FontSize',18)

