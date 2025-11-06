% Analysis of all vertical profiles from non-precipitating clouds measured
% during the VOCALS-REx campaign on-board the NCAR C130 aircract

% By Andrew John Buggee


%% Plot the ensemble distribution of effective radius and optical depth
% for non-precipitating clouds. 
% ----- For Vertical Profiles -----


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true



elseif strcmp(whatComputer, 'andrewbuggee')==true

    % Location of ensemble data
    folderpath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/VOCALS_REx/',...
        'vocals_rex_data/NCAR_C130/SPS_1/'];



    % --- non-precip profiles only, LWC>0.03, Nc>25 ----
    % load([folderpath,...
    % 'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_30-Oct-2025_rev1.mat'])

    load([folderpath,...
      'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_drizzleLWP-threshold_5_02-Nov-2025.mat'])


end





%%
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth to normalize each
%  profile

% In order to compute a mean vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. 

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);

% Store all cloud optical depths
tau_c = zeros(1, length(ensemble_profiles));


normalized_altitude = cell(1, length(ensemble_profiles));




for nn = 1:length(ensemble_profiles)

    % first we need to normalize the vertical component of all profiles
    normalized_altitude{nn} = (ensemble_profiles{nn}.altitude - min(ensemble_profiles{nn}.altitude))./...
        (max(ensemble_profiles{nn}.altitude) - min(ensemble_profiles{nn}.altitude));

    % the data is stored in altitude space.

    if ensemble_profiles{nn}.flag_2DC_data_is_conforming==true

        re = ensemble_profiles{nn}.re;
    else

        % Since we are analyzing non-precipitating clouds, as long as the
        % 2DC LWP is small, we can simply take the CDP data as the
        % effective radius data
        if ensemble_profiles{nn}.lwp_2DC<5

            re = ensemble_profiles{nn}.re_CDP;

        end

    end

    lwc = ensemble_profiles{nn}.lwc;
    Nc_total = sum(ensemble_profiles{nn}.Nc, 1);

    % check to see if any values are 0
    if sum(re==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        re(re==0) = 0.01; % microns;

    end

    if sum(lwc==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        lwc(lwc==0) = 0.001; % microns;

    end

    if sum(Nc_total==0)>=1

        % This happens sometimes in multi layer clouds. Set it to be a
        % small non-zero value
        Nc_total(Nc_total==0) = 0.1; % microns;

    end



    % Also store the cloud optical thickness for each profile
    tau_c(nn) = max(ensemble_profiles{nn}.tau);

    % Check if there is ever a 0
    if tau_c(nn)==0

        error('Why?')

    end




    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_altitude{nn}>=bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_altitude{nn}>bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; reshape(re(index_segment),[],1)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; reshape(lwc(index_segment), [],1)];

        % store the total number concentration
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; reshape(Nc_total(index_segment), [],1)];



    end



end



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% We know want to determine the statistics of the effective radius, lwc and
% total number concentration as a function of height within the cloud. At
% each normalized level, let's fit a normal, lognormal and gamma function
% to the distribution and see which is the best fit by check the
% chi-squared test

significance_lvl = 0.05;

% store the refection of each null hypothesis and the p-value for each
% chi-squared test

re_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
re_p_normal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
re_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
re_p_gamma = zeros(1, size(vertically_segmented_attributes,1));


lwc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));



Nc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));

for bb = 1:size(vertically_segmented_attributes, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    re_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'normal');
    
    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_normal(bb), re_p_normal(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb),...
        'alpha', significance_lvl);

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'lognormal');

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb),...
        'alpha', significance_lvl);

    % fit the effective radius data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    re_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(vertically_segmented_attributes{bb,1});

    % *** Chi-squared doesn't work well ***
    % chi-squared method to fail. Try using the Kolmogorov-Smirnov test,
    % which is more robust to outliers
    % [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
    %     'alpha', significance_lvl, 'NParams', 2);
    [re_reject_gamma(bb), re_p_gamma(bb)] = kstest(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb),...
        'alpha', significance_lvl);

    % make plot of normal log-normal and gamma fits
        % figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
        % subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
        % subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
        % set(gcf, 'Position', [0 0 1200 500])








    % -------------------------------------------
    % ------- LIQUID WATER CONTENT FITTING ------
    % -------------------------------------------


    % fit the liquid water content data to a normal distribution
    lwc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'normal');
    [lwc_reject_normal(bb), lwc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_normal(bb),...
        'alpha', significance_lvl, 'NParams', 2);


    % fit the liquid water content data to a log-normal distribution
    lwc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'lognormal');
    [lwc_reject_lognormal(bb), lwc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_lognormal(bb),...
        'alpha', significance_lvl, 'NParams', 2);

    % fit the liquid water content data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    lwc_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(vertically_segmented_attributes{bb,2});
    [lwc_reject_gamma(bb), lwc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF', lwc_fit_gamma(bb),...
        'alpha', significance_lvl, 'NParams', 2);

    % make plot of normal log-normal and gamma fits
    %     figure; subplot(1,3,1); plot(lwc_fit_normal(bb)); title('Normal Fit'); xlabel('LWC (g/m^{2})')
    %     subplot(1,3,2); plot(lwc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (g/m^{2})')
    %     subplot(1,3,3); plot(lwc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (g/m^{2})')
    %     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- NUMBER CONCENTRATION FITTING ------
    % -------------------------------------------


    % fit the number concentration data to a normal distribution
    Nc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'normal');
    [Nc_reject_normal(bb), Nc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_normal(bb),...
        'alpha', significance_lvl, 'NParams', 2);

    % fit the number concentration content data to a log-normal distribution
    Nc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'lognormal');
    [Nc_reject_lognormal(bb), Nc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_lognormal(bb),...
        'alpha', significance_lvl, 'NParams', 2);

    % fit the total number concentration data to a gamma distribution - use my custom
    % libRadtran gamma distribution
    Nc_fit_gamma(bb) = prob.GammaDistribution_libRadtran.fit(vertically_segmented_attributes{bb,3});
    [Nc_reject_gamma(bb), Nc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF', Nc_fit_gamma(bb),...
        'alpha', significance_lvl, 'NParams', 2);

    % make plot of normal log-normal and gamma fits
    %     figure; subplot(1,3,1); plot(Nc_fit_normal(bb)); title('Normal Fit'); xlabel('N_c (cm^{-3})')
    %     subplot(1,3,2); plot(Nc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (cm^{-3})')
    %     subplot(1,3,3); plot(Nc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (cm^{-3})')
    %     set(gcf, 'Position', [0 0 1200 500])


end



% -------------------------------------------
% ----------- OPTICAL DEPTH FITTING ---------
% -------------------------------------------
% Let's also fit these three distributions to the optical depth data

% fit the number concentration data to a normal distribution
tauC_fit_normal = fitdist(tau_c, 'normal');
[tauC_reject_normal, tauC_p_normal] = chi2gof(tau_c, 'CDF', tauC_fit_normal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the number concentration content data to a log-normal distribution
tauC_fit_lognormal = fitdist(tau_c, 'lognormal');
[tauC_reject_lognormal, tauC_p_lognormal] = chi2gof(tau_c, 'CDF', tauC_fit_lognormal,...
    'alpha', significance_lvl, 'NParams', 2);

% fit the total number concentration data to a gamma distribution - use my custom
% libRadtran gamma distribution
tauC_fit_gamma = prob.GammaDistribution_libRadtran.fit(tau_c);
[tauC_reject_gamma, tauC_p_gamma] = chi2gof(tau_c, 'CDF', tauC_fit_gamma,...
    'alpha', significance_lvl, 'NParams', 2);






% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
[max_re_p, idx_re_p] = max([re_p_normal; re_p_lognormal; re_p_gamma],[], 1, "omitnan");

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]);
% title('r_e best distribution fit'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1, "omitnan");

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]);
% title('LWC best distribution fit'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1, "omitnan");

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]);
% title('N_c best distribution fit'); ylabel('Counts')



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the mean LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_mean = zeros(n_bins, 1);


% ---- most common best fit distribution for LWC was is the normal dist ---
lwc_mean = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);


% ---- most common best fit distribution for N_c was is the normal dist ---
Nc_mean = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);




for bb = 1:n_bins


    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_mean(bb) = exp(re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2 /2);

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));



    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % compute the mean value for the current bin
    % the mean of the distribution (the standard way of computing the expected value)
    % is also the mean of the normal distribution. They are identical.
    lwc_mean(bb) = mean(vertically_segmented_attributes{bb,2});       % g/cm^3 - mean liqiud water content
    %lwc_mean(bb) = lwc_fit_gamma(bb).mean;       % g/cm^3 - mean liqiud water content

    % compute the standard deviation of the current bin
    % the std of the distribution (the standard way of computing the squareroot of the variance)
    % is also the std of the normal distribution. They are identical.
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % g/cm^3 - standard deviation
    %lwc_std(bb) = lwc_fit_gamma(bb).std;         % g/cm^3 - standard deviation



    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

    % compute the mean value for the current bin
    %Nc_mean(bb) = Nc_fit_gamma(bb).mean;       % cm^(-3) - mean number concentration
    Nc_mean(bb) = Nc_fit_normal(bb).mean;       % cm^(-3) - mean number concentration

    % compute the standard deviation of the current bin
    %Nc_std(bb) = Nc_fit_gamma(bb).std;         % cm^(-3) - standard deviation
    Nc_std(bb) = Nc_fit_normal(bb).std;         % cm^(-3) - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% -------------- Make a subplot of all 3 mean profiles -----------------
% ----------------------------------------------------------------------




figure;

% plot the mean effective radius
subplot(1,3,1)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_logNormal_mean - re_logNormal_std; flipud(re_logNormal_mean + re_logNormal_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(re_logNormal_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex', 'Fontsize',30)
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize',30)

% set x axis boundaries
xlim([4, 10])                   % microns



% plot the mean liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex', 'FontSize', 30)
ylabel('Normalized Altitude', 'Interpreter','latex', 'FontSize', 30)

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

% set the figure title
title(['Mean Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')





% plot the mean droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_mean-Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(Nc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 700])




% ----------------------------------------------------------------------
% ----------------------- Adiabatic Curve Fits -------------------------
% ----------------------------------------------------------------------

% ------------------------- LIQUID WATER CONTENT ---------------------------
% ----- Fit an adiabatic curve to the mean liquid water content profile -----
nudge_from_top =2;
lwc_slope = (lwc_mean(end-nudge_from_top) - lwc_mean(1))/(bin_center(end-nudge_from_top) - bin_center(1));
lwc_intercept = lwc_mean(1) - lwc_slope*bin_center(1);
lwc_adiabatic_fit = lwc_slope*bin_center + lwc_intercept;

% add to subplot(1,3,1)
subplot(1,3,2); hold on
plot(lwc_adiabatic_fit(1:end-nudge_from_top), bin_center(1:end-nudge_from_top), 'k', 'LineWidth',2)

% -- Include a legend in the lower right-hand corner of the 2nd subplot --
legend({'Standard Deviation', 'Mean Profile', 'Adiabatic Fit'}, 'Interpreter','latex',...
    'Location','southeast', 'FontSize', 17)



% -------------------- EFFECTIVE DROPLET RADIUS -----------------------
% Plot an adiabatic curve fit to the mean droplet radius profile
nudge_from_bottom = 3;

% ----- Fit an adiabatic curve to the mean droplet profile -----
% use droplet profile function to create adiabatic fit
re_adiabatic_fit = create_droplet_profile2([re_logNormal_mean(end), re_logNormal_mean(1 + nudge_from_bottom)],...
    bin_center(1+nudge_from_bottom:end),'altitude', 'adiabatic');

% derive droplet profile from adiabatic fit of LWC


% add to subplot(1,3,1)
subplot(1,3,1); hold on
plot(re_adiabatic_fit, bin_center(1+nudge_from_bottom:end), 'k', 'LineWidth',2)





% ---------------------- DROPLET NUMBER CONCENTRATION ---------------------------
% ----- Fit an adiabatic curve to Number Concentration based on Adiabatic fits above -----
% Nc_adiabatic_fit = 3*(lwc_adiabatic_fit./1e6)./...
%                 (4*pi*(re_adiabatic_fit*1e-4).^3);




% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<1 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');



