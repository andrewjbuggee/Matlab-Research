%% TITLE: EPCAPE_Pier_CCNclosure_and_correction.m
% 
% PURPOSE: Performs bottom-up CCN closure using composition, size 
% distribution, and CCN measurements.
%
% INPUTS: EPCAPE data (.nc)
% 
% OUTPUTS: .MAT files 
% 
% AUTHOR: Jeramy Dedrick
%         Scripps Institution of Oceanography, La Jolla, CA
%         16 April 2025
%         Edited: 3 November 2025

format long
clear all; close all; clc


%% Load data

load('EPCAPE_Pier_kappa4CCN.mat')

matlab_colors = [0, 0.4470, 0.7410;
          	    0.8500, 0.3250, 0.0980;        	
          	    0.9290, 0.6940, 0.1250;
          	    0.4940, 0.1840, 0.5560;        	
          	    0.4660, 0.6740, 0.1880;
          	    0.3010, 0.7450, 0.9330;          	
          	    0.6350, 0.0780, 0.1840];



%% Calculate Critical Diameter from using Kappa-Kohler Theory (P&K 2006)
close all;clc

% thermodynamic constants
sigma_sa = 0.072; % surface tension of the solution/air interface, J/m2
M_w      = 0.018015; % molecular weight of water, kg/mol
rho_w    = 997.1; % density of water, kg/m3
R        = 8.315; % universal gas constant, J/mol/K
T        = 298.15; % temperature

% CCN supersaturations
CCN_supersat = [0.1 0.2 0.4 0.6 0.8 1.0];

Dd_ACSM = logspace(log10(25),log10(300),100) .* (1e-9);   % dry diameter within ACSM measurement range, m
D_wet   = logspace(log10(10),log10(50000),100) .* (1e-9); % wet diameter that covers relevant size range, m

% empty matrices to fill
Pier_Dc_kappa_SSDust_Berta     = NaN(length(CCN_supersat), length(EPCAPE_time_30min));
Pier_Dc_kappa_noSSnoDust_Berta = NaN(length(CCN_supersat), length(EPCAPE_time_30min));


% perform loop to compute critical diameters using kappa (with/without sea
% salt and dust)
for i = 1:length(EPCAPE_time_30min)

    % empty matrices to fill
    S_d_kappa_SSDust_Berta       = NaN(length(D_wet), length(Dd_ACSM));
    S_d_kappa_noSSnoDust_Berta   = NaN(length(D_wet), length(Dd_ACSM));
    
    Smax_kappa_SSDust_Berta      = NaN(1, length(D_wet)); 
    Smax_kappa_noSSnoDust_Berta  = NaN(1, length(D_wet)); 


    for j = 1:length(Dd_ACSM)
     
        for k = 1:length(D_wet)
    
            if D_wet(k) < Dd_ACSM(j) % if wet diameter smaller than dry

                S_d_kappa_SSDust_Berta(k,j)      = NaN;
                S_d_kappa_noSSnoDust_Berta(k,j) = NaN;
    
            else
               % EQ. (6) Petters and Kreidenweis (2006)
               S_d_kappa_SSDust_Berta(k,j) = ((D_wet(k).^3 - Dd_ACSM(j).^3) ./ ...
                                             (D_wet(k).^3 - (Dd_ACSM(j).^3 .* (1 - Berta_pier_kappa_SSDust_retime_30min(i))))) .* ...
                                             exp((4 .* sigma_sa .* M_w) ./ (R .* T .* rho_w .* D_wet(k)));

               S_d_kappa_noSSnoDust_Berta(k,j) = ((D_wet(k).^3 - Dd_ACSM(j).^3) ./ ...
                                                  (D_wet(k).^3 - (Dd_ACSM(j).^3 .* (1 - Berta_pier_kappa_noSSnoDust_retime_30min(i))))) .* ...
                                                  exp((4 .* sigma_sa .* M_w) ./ (R .* T .* rho_w .* D_wet(k)));
            end
    
        end  
    
        % find index of "peak" (maximum) saturation ratio for each dry diameter
        get_max_S_kappa_SSDust_Berta     = find(S_d_kappa_SSDust_Berta(:,j) == max(S_d_kappa_SSDust_Berta(:,j)));
        get_max_S_kappa_noSSnoDust_Berta = find(S_d_kappa_noSSnoDust_Berta(:,j) == max(S_d_kappa_noSSnoDust_Berta(:,j)));

        % get peak saturation ratio and convert to percent
        if isempty(get_max_S_kappa_SSDust_Berta)
            Smax_kappa_SSDust_Berta(j) = NaN;
        else
            Smax_kappa_SSDust_Berta(j) = (S_d_kappa_SSDust_Berta(get_max_S_kappa_SSDust_Berta, j) - 1) .* 100;
        end

        if isempty(get_max_S_kappa_noSSnoDust_Berta)
            Smax_kappa_noSSnoDust_Berta(j) = NaN;
        else
            Smax_kappa_noSSnoDust_Berta(j) = (S_d_kappa_noSSnoDust_Berta(get_max_S_kappa_noSSnoDust_Berta, j) - 1) .* 100;
        end


    end
    

    % interpolate to find the Dd that pertains to Smax=SCCN
    if all(isnan(Smax_kappa_SSDust_Berta))
        Pier_Dc_kappa_SSDust_Berta(:,i) = NaN(length(CCN_supersat),1);
    else
        Pier_Dc_kappa_SSDust_Berta(:,i) = interp1(Smax_kappa_SSDust_Berta, ...
                                                Dd_ACSM, ...
                                                CCN_supersat') .* (1e6);
    end

    if all(isnan(Smax_kappa_noSSnoDust_Berta))
        Pier_Dc_kappa_noSSnoDust_Berta(:,i) = NaN(length(CCN_supersat),1);
    else
        Pier_Dc_kappa_noSSnoDust_Berta(:,i) = interp1(Smax_kappa_noSSnoDust_Berta, ...
                                                          Dd_ACSM, ...
                                                          CCN_supersat') .* (1e6);
    end


    disp(strcat('calculating Dc: ', num2str(i), '/', num2str(length(EPCAPE_time_30min))))

end



%% Predicted CCN from Retrieved Dcrit
close all; clc

Pier_NCCN_pred_kappa_SSDust_Berta      = NaN(length(CCN_supersat), length(EPCAPE_time_30min));
Pier_NCCN_pred_kappa_noSSnoDust_Berta  = NaN(length(CCN_supersat), length(EPCAPE_time_30min));

for i = 1:length(EPCAPE_time_30min)

    for j = 1:length(CCN_supersat)

        % find SMPS diameter closest to and greater than critical diameter
        find_Dc_kappa_SSDust_Berta     = find(SMPS_d_mid >= Pier_Dc_kappa_SSDust_Berta(j,i), 1, 'first');
        find_Dc_kappa_noSSnoDust_Berta = find(SMPS_d_mid >= Pier_Dc_kappa_noSSnoDust_Berta(j,i), 1, 'first');

        % sum over number concentration to get CCN 
        if isempty(find_Dc_kappa_SSDust_Berta)
            Pier_NCCN_pred_kappa_SSDust_Berta(j,i) = NaN;
        else
            Pier_NCCN_pred_kappa_SSDust_Berta(j,i) = nansum(SMPS_NSD_retime_30min(find_Dc_kappa_SSDust_Berta:end,i));            
        end

        if isempty(find_Dc_kappa_noSSnoDust_Berta)
            Pier_NCCN_pred_kappa_noSSnoDust_Berta(j,i) = NaN;
        else
            Pier_NCCN_pred_kappa_noSSnoDust_Berta(j,i) = nansum(SMPS_NSD_retime_30min(find_Dc_kappa_noSSnoDust_Berta:end,i));            
        end

    end
end

% omit when predicted value is zero
Pier_NCCN_pred_kappa_SSDust_Berta(Pier_NCCN_pred_kappa_SSDust_Berta==0)         = NaN;
Pier_NCCN_pred_kappa_noSSnoDust_Berta(Pier_NCCN_pred_kappa_noSSnoDust_Berta==0) = NaN;



%% Calculate ratio of measured to predicted

% measured over predicted ratio
ratio_pred_kappa_SSDust_Berta = [CCNcolaavg_spectra_retime_30min(1,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(1,:);
                                 CCNcolaavg_spectra_retime_30min(2,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(2,:);
                                 CCNcolaavg_spectra_retime_30min(3,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(3,:);
                                 CCNcolaavg_spectra_retime_30min(4,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(4,:);
                                 CCNcolaavg_spectra_retime_30min(5,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(5,:);
                                 CCNcolaavg_spectra_retime_30min(6,:) ./ Pier_NCCN_pred_kappa_SSDust_Berta(6,:)];
inverse_ratio_pred_kappa_noSSnoDust_Berta = 1./ratio_pred_kappa_SSDust_Berta;

ratio_pred_kappa_noSSnoDust_Berta = [CCNcolaavg_spectra_retime_30min(1,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(1,:);
                                     CCNcolaavg_spectra_retime_30min(2,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(2,:);
                                     CCNcolaavg_spectra_retime_30min(3,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(3,:);
                                     CCNcolaavg_spectra_retime_30min(4,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(4,:);
                                     CCNcolaavg_spectra_retime_30min(5,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(5,:);
                                     CCNcolaavg_spectra_retime_30min(6,:) ./ Pier_NCCN_pred_kappa_noSSnoDust_Berta(6,:)];


%% Bad Data based on Ratio being too high

% threshold for unrealistic ratio
bad_thresh = 4;

for i = 1:6

    remove_SSDust_idx     = [];
    remove_noSSnoDust_idx = [];
    remove_SSDust_idx2    = [];

    remove_SSDust_idx     = find(ratio_pred_kappa_SSDust_Berta(i,:) > bad_thresh);
    remove_SSDust_idx2    = find(inverse_ratio_pred_kappa_noSSnoDust_Berta(i,:) > bad_thresh);
    remove_noSSnoDust_idx = find(ratio_pred_kappa_noSSnoDust_Berta(i,:) > bad_thresh);

    Pier_NCCN_pred_kappa_SSDust_Berta(i,remove_SSDust_idx) = NaN;
    ratio_pred_kappa_SSDust_Berta(i,remove_SSDust_idx) = NaN;
    inverse_ratio_pred_kappa_noSSnoDust_Berta(i,remove_SSDust_idx2) = NaN;

    Pier_NCCN_pred_kappa_noSSnoDust_Berta(i,remove_noSSnoDust_idx) = NaN;
    ratio_pred_kappa_noSSnoDust_Berta(i,remove_noSSnoDust_idx) = NaN;

end


%% Some Time Variables
close all; clc

EPCAPE_time_30min_floor   = floor(datenum(EPCAPE_time_30min));
EPCAPE_time_30min_datevec = datevec(EPCAPE_time_30min);
dayhours                  = 0:23;

% create time grids
startTime            = datetime(2023, 2, 15, 0, 0, 0);
endTime              = datetime(2024, 2, 15, 0, 0, 0);
EPCAPE_time_1hr      = startTime:hours(1):endTime;
EPCAPE_time_3hr      = startTime:hours(3):endTime;
EPCAPE_time_1day     = startTime:days(1):endTime;
EPCAPE_time_7day     = startTime:days(7):endTime;


%% Find Average Ratio for a Week 
close all; clc

EPCAPE_weeks         = EPCAPE_time_30min(1):days(7):EPCAPE_time_30min(end);
EPCAPE_weeks_datenum = datenum(EPCAPE_weeks)';

% empty matrices
ratio_pred_kappa_SSDust_Berta_weekly        = NaN(6,length(EPCAPE_weeks),2);
ratio_pred_kappa_noSSnoDust_Berta_weekly    = NaN(6,length(EPCAPE_weeks),2);
f_org_noSSnoDust_weekly                     = NaN(1,length(EPCAPE_weeks));
f_NH4NO3_noSSnoDust_weekly                  = NaN(1,length(EPCAPE_weeks));
f_NH42SO4_noSSnoDust_weekly                 = NaN(1,length(EPCAPE_weeks));
f_NH4HSO4_noSSnoDust_weekly                 = NaN(1,length(EPCAPE_weeks));
f_org_SSDust_weekly                         = NaN(1,length(EPCAPE_weeks));
f_NH4NO3_SSDust_weekly                      = NaN(1,length(EPCAPE_weeks));
f_NH42SO4_SSDust_weekly                     = NaN(1,length(EPCAPE_weeks));
f_NH4HSO4_SSDust_weekly                     = NaN(1,length(EPCAPE_weeks));
f_seasalt_SSDust_weekly                     = NaN(1,length(EPCAPE_weeks));
f_dust_SSDust_weekly                        = NaN(1,length(EPCAPE_weeks));
ratio_pred_kappa_SSDust_Berta_2hrly_weekly  = NaN(6,12,length(EPCAPE_weeks));
ratio_pred_kappa_SSDust_Berta_2hrly_weekly_time = NaN(12,length(EPCAPE_weeks));

for i = 1:length(EPCAPE_weeks)-1
 
    if i == 32 % skip this week (no data)

    else
    
    % get indices for the specific week
    get_week = find(EPCAPE_time_30min_floor >= EPCAPE_weeks_datenum(i) & ...
                    EPCAPE_time_30min_floor < EPCAPE_weeks_datenum(i+1));
    get_week_datevec = EPCAPE_time_30min_datevec(get_week,:);

    % get indices for a random day in a week
    get_random_day = randi([min(get_week_datevec(:,3)) max(get_week_datevec(:,3))]);
    get_day        = find(get_week_datevec(:,3)==get_random_day);

  
    % get data for the specific weeek
    ratio_pred_kappa_SSDust_Berta_get_week      = ratio_pred_kappa_SSDust_Berta(:,get_week);
    ratio_pred_kappa_noSSnoDust_Berta_get_week  = ratio_pred_kappa_noSSnoDust_Berta(:,get_week);
    f_org_noSSnoDust_get_week                   = Berta_pier_f_org_noSSnoDust_retime_30min(:,get_week);
    f_NH4NO3_noSSnoDust_get_week                = Berta_pier_f_NH4NO3_noSSnoDust_retime_30min(:,get_week);
    f_NH42SO4_noSSnoDust_get_week               = Berta_pier_f_NH42SO4_noSSnoDust_retime_30min(:,get_week);
    f_NH4HSO4_noSSnoDust_get_week               = Berta_pier_f_NH4HSO4_noSSnoDust_retime_30min(:,get_week);    
    f_org_SSDust_get_week                       = Berta_pier_f_org_SSDust_retime_30min(:,get_week);
    f_NH4NO3_SSDust_get_week                    = Berta_pier_f_NH4NO3_SSDust_retime_30min(:,get_week);
    f_NH42SO4_SSDust_get_week                   = Berta_pier_f_NH42SO4_SSDust_retime_30min(:,get_week);
    f_NH4HSO4_SSDust_get_week                   = Berta_pier_f_NH4HSO4_SSDust_retime_30min(:,get_week);
    f_seasalt_SSDust_get_week                   = Berta_pier_f_seasalt_SSDust_retime_30min(:,get_week);
    f_dust_SSDust_get_week                      = Berta_pier_f_dust_SSDust_retime_30min(:,get_week);

    % get statistics the specific week
    ratio_pred_kappa_SSDust_Berta_weekly(:,i,1)     = nanmean(ratio_pred_kappa_SSDust_Berta_get_week,2);
    ratio_pred_kappa_noSSnoDust_Berta_weekly(:,i,1) = nanmean(ratio_pred_kappa_noSSnoDust_Berta_get_week,2);
    ratio_pred_kappa_SSDust_Berta_weekly(:,i,2)     = nanstd(ratio_pred_kappa_SSDust_Berta_get_week,0,2);
    ratio_pred_kappa_noSSnoDust_Berta_weekly(:,i,2) = nanstd(ratio_pred_kappa_noSSnoDust_Berta_get_week,0,2);
    f_org_noSSnoDust_weekly(:,i)                    = nanmean(Berta_pier_f_org_noSSnoDust_retime_30min(:,get_week));
    f_NH4NO3_noSSnoDust_weekly(:,i)                 = nanmean(Berta_pier_f_NH4NO3_noSSnoDust_retime_30min(:,get_week));
    f_NH42SO4_noSSnoDust_weekly(:,i)                = nanmean(Berta_pier_f_NH42SO4_noSSnoDust_retime_30min(:,get_week));
    f_NH4HSO4_noSSnoDust_weekly(:,i)                = nanmean(Berta_pier_f_NH4HSO4_noSSnoDust_retime_30min(:,get_week));    
    f_org_SSDust_weekly(:,i)                        = nanmean(Berta_pier_f_org_SSDust_retime_30min(:,get_week));
    f_NH4NO3_SSDust_weekly(:,i)                     = nanmean(Berta_pier_f_NH4NO3_SSDust_retime_30min(:,get_week));
    f_NH42SO4_SSDust_weekly(:,i)                    = nanmean(Berta_pier_f_NH42SO4_SSDust_retime_30min(:,get_week));
    f_NH4HSO4_SSDust_weekly(:,i)                    = nanmean(Berta_pier_f_NH4HSO4_SSDust_retime_30min(:,get_week));
    f_seasalt_SSDust_weekly(:,i)                    = nanmean(Berta_pier_f_seasalt_SSDust_retime_30min(:,get_week));
    f_dust_SSDust_weekly(:,i)                       = nanmean(Berta_pier_f_dust_SSDust_retime_30min(:,get_week));


    % get data for the random day in the week
    ratio_day_timetable = timetable(datetime(get_week_datevec(get_day,:)), ...
                                     ratio_pred_kappa_noSSnoDust_Berta_get_week(:,get_day)');
    if isempty(ratio_day_timetable)
    else
    ratio_day_2hrly = retime(ratio_day_timetable, 'regular', 'mean', 'TimeStep', hours(2));
    ratio_pred_kappa_SSDust_Berta_2hrly_weekly(:,:,i) = ratio_day_2hrly.Var1';
    ratio_pred_kappa_SSDust_Berta_2hrly_weekly_time(:,i) = datenum(ratio_day_2hrly.Time);
    end
    
    end

end


%% Plot Weekly Ratio - noSSnoDust


f_vol_bar = [f_org_noSSnoDust_weekly', ...
             f_NH42SO4_noSSnoDust_weekly', ...
             f_NH4HSO4_noSSnoDust_weekly',...
             f_NH4NO3_noSSnoDust_weekly'];

f_vol_color = [51 160 44;
               227 26 28;
               171 77 78;
               31 120 180] ./ 255;

% f_vol_barchart = bar(f_vol_bar, 'stacked');

figure(1)
set(gcf, 'Position', [100 100 750 400])
subplot('Position', [0.08 0.12 0.7 0.83])
for i = 4:6
    for j = 14:size(ratio_pred_kappa_noSSnoDust_Berta_weekly, 2)
   
        colororder({'k', 'k'})

        yyaxis right
        f_vol_barchart = bar(j, f_vol_bar(j,:), 'stacked');
        for k = 1:length(f_vol_barchart)
            f_vol_barchart(k).FaceColor = f_vol_color(k,:);
            f_vol_barchart(k).FaceAlpha = 1.0;
        end
        hold on
        ylim([0 3])
        yticks([0:0.2:1])
        xticks([14:4:50])
        xlim([13 53])
        


        yyaxis left
        errorbar(j, ...
                 ratio_pred_kappa_noSSnoDust_Berta_weekly(i,j,1), ...
                 ratio_pred_kappa_noSSnoDust_Berta_weekly(i,j,2), 'Color', matlab_colors(i,:))
        hold on
        plot(j, ratio_pred_kappa_noSSnoDust_Berta_weekly(i,j,1), '.', 'Color',  matlab_colors(i,:), 'MarkerSize',20)
        
    end
end
yline(0)
hold off
text(56, -0.4, 'f_{vol}', 'Rotation', 90, 'FontSize', 10)
text(58, -0.1, 'NR-Org', 'Color', f_vol_color(1,:), 'FontSize',8)
text(58, -0.2, 'NR-(NH_{4})_{2}SO_{4}', 'Color', f_vol_color(2,:), 'FontSize',8)
text(58, -0.3, 'NR-(NH_{4})HSO_{4}', 'Color', f_vol_color(3,:), 'FontSize',8)
text(58, -0.4, 'NR-(NH_{4})NO_{3}', 'Color', f_vol_color(4,:), 'FontSize',8)
text(54, 1, '0.6% Ratio', 'Color', matlab_colors(4,:), 'FontSize', 8)
text(54, 0.9, '0.8% Ratio', 'Color', matlab_colors(5,:), 'FontSize', 8)
text(54, 0.8, '1.0% Ratio', 'Color', matlab_colors(6,:), 'FontSize', 8)
xlabel('Week #')
ylabel('CCN_{meas.}/CCN_{pred.}')
ylim([-0.6 1.2])
yticks([0:0.2:1.2])
xticks([14:4:50])
xlim([13 53])
set(gca, 'FontSize', 10)


%% Plot Weekly Ratio - SSDust
clc

f_vol_bar = [f_org_SSDust_weekly', ...
             f_NH42SO4_SSDust_weekly', ...
             f_NH4HSO4_SSDust_weekly',...
             f_NH4NO3_SSDust_weekly',...
             f_seasalt_SSDust_weekly',...
             f_dust_SSDust_weekly'];


f_vol_color = [51 160 44;
               227 26 28;
               171 77 78;
               31 120 180;
               235 89 186;
               150 75 0] ./ 255;

% f_vol_barchart = bar(f_vol_bar, 'stacked');

figure(2)
set(gcf, 'Position', [100 100 750 400])
subplot('Position', [0.08 0.12 0.7 0.83])
for i = 4:6
    for j = 14:size(ratio_pred_kappa_SSDust_Berta_weekly, 2)
   
        colororder({'k', 'k'})

        yyaxis right
        f_vol_barchart = bar(j, f_vol_bar(j,:), 'stacked');
        for k = 1:length(f_vol_barchart)
            f_vol_barchart(k).FaceColor = f_vol_color(k,:);
            f_vol_barchart(k).FaceAlpha = 1.0;
        end
        hold on
        ylim([0 3])
        yticks([0:0.2:1])
        xticks([14:4:50])
        xlim([13 53])
        


        yyaxis left
        errorbar(j, ...
                 ratio_pred_kappa_SSDust_Berta_weekly(i,j,1), ...
                 ratio_pred_kappa_SSDust_Berta_weekly(i,j,2), 'Color', matlab_colors(i,:))
        hold on
        plot(j, ratio_pred_kappa_SSDust_Berta_weekly(i,j,1), '.', 'Color',  matlab_colors(i,:), 'MarkerSize',20)
        
    end
end
yline(0)
hold off
text(56, -0.4, 'f_{vol}', 'Rotation', 90, 'FontSize', 10)
text(58, -0.1, 'NR-Org', 'Color', f_vol_color(1,:), 'FontSize',8)
text(58, -0.2, 'NR-(NH_{4})_{2}SO_{4}', 'Color', f_vol_color(2,:), 'FontSize',8)
text(58, -0.3, 'NR-(NH_{4})HSO_{4}', 'Color', f_vol_color(3,:), 'FontSize',8)
text(58, -0.4, 'NR-(NH_{4})NO_{3}', 'Color', f_vol_color(4,:), 'FontSize',8)
text(58, -0.5, 'sea salt', 'Color', f_vol_color(5,:), 'FontSize',8)
text(58, -0.6, 'dust', 'Color', f_vol_color(6,:), 'FontSize',8)
text(54, 1, '0.6% Ratio', 'Color', matlab_colors(4,:), 'FontSize', 8)
text(54, 0.9, '0.8% Ratio', 'Color', matlab_colors(5,:), 'FontSize', 8)
text(54, 0.8, '1.0% Ratio', 'Color', matlab_colors(6,:), 'FontSize', 8)
xlabel('Week #')
ylabel('CCN_{meas.}/CCN_{pred.}')
ylim([-0.6 1.2])
yticks([0:0.2:1.2])
xticks([14:4:50])
xlim([13 53])
set(gca, 'FontSize', 10)



%% Calculate Campaign Scaling Factors
clc

low_salt_weeks = 14:31;
hi_salt_weeks  = 32:53;

ScaleFactor_kappa_SSDust_campaign     = [1 ./ nanmean(ratio_pred_kappa_SSDust_Berta_weekly(4:6,:,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_SSDust_Berta_weekly(4:6,:,1), 0,'all')];
ScaleFactor_kappa_noSSnoDust_campaign = [1 ./ nanmean(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,:,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,:,1), 0,'all')];

ScaleFactor_kappa_SSDust_low_salt     = [1 ./ nanmean(ratio_pred_kappa_SSDust_Berta_weekly(4:6,low_salt_weeks,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_SSDust_Berta_weekly(4:6,low_salt_weeks,1), 0,'all')];
ScaleFactor_kappa_noSSnoDust_low_salt = [1 ./ nanmean(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,low_salt_weeks,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,low_salt_weeks,1), 0,'all')];

ScaleFactor_kappa_SSDust_hi_salt     = [1 ./ nanmean(ratio_pred_kappa_SSDust_Berta_weekly(4:6,hi_salt_weeks,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_SSDust_Berta_weekly(4:6,hi_salt_weeks,1), 0,'all')];
ScaleFactor_kappa_noSSnoDust_hi_salt = [1 ./ nanmean(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,hi_salt_weeks,1), 'all'), ...
                                         nanstd(ratio_pred_kappa_noSSnoDust_Berta_weekly(4:6,hi_salt_weeks,1), 0,'all')];


% correct CCN measurements by low_salt scaling factor
Pier_CCN_corrected_kappa_SSDust_low_salt_30min = CCNcolaavg_spectra_retime_30min .* ScaleFactor_kappa_noSSnoDust_low_salt(1,1);

