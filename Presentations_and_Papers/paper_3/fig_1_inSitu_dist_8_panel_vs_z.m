%% Combine VOCALS-REx and ORACLES - 25-75 filled band, 10th/90th as dashed lines
%
% The cell blocks above display "spread" in two ways that both hide the
% asymmetry of the underlying lognormal-like distributions:
%
%   1) IQR version:         x = [median - IQR/2;  flipud(median + IQR/2)]
%   2) Average-deviation:   symmetric by construction if the two one-sided
%                           deviations happen to be equal; and more subtly,
%                           both approaches center the band ON the median
%                           rather than using percentile values directly.
%
% The IQR case is the clearest example: the 25th-75th percentile range is
% split evenly on either side of the median regardless of where Q1 and Q3
% actually lie, so the resulting band is symmetric around the median by
% construction. Skewness in the distribution is thrown away.
%
% Here we plot percentiles directly:
%   filled band  : 25th to 75th percentile
%   dashed lines : 10th and 90th percentiles (same color as the median)
%   solid line   : median (50th percentile)
%
% Because each edge of the band and each dashed line is the actual
% percentile of the sampled data at that vertical level, the left/right
% widths relative to the median reflect the true asymmetry in the
% distribution.


% -------------------------------------------------------------------------
% ------ VOCALS-REx: load, separate drizzle, compute percentiles ----------
% -------------------------------------------------------------------------
clear variables

tauC_limit = 3;

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

end

load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

n_bins     = 30;
bin_edges  = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

vertSeg_nd = cell(n_bins, 4);
vertSeg_d  = cell(n_bins, 4);

num_removed = 0;
num_kept    = 0;

for nn = 1:N_profiles

    if max(ensemble_profiles{nn}.tau) < tauC_limit
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end

    % cloud profiles should range from cloud top to cloud bottom
    z_norm   = (ensemble_profiles{nn}.altitude - min(ensemble_profiles{nn}.altitude)) ./...
        (max(ensemble_profiles{nn}.altitude) - min(ensemble_profiles{nn}.altitude));

    significance_lvl = 0.1;
    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nc, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);
    v_eff_bin = 1 ./ (3 + gammaFit.alpha);

    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re);
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re_CDP);
        end
        lwc_bin   = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin    = fliplr(ensemble_profiles{nn}.total_Nc);
        v_eff_bin = flipud(v_eff_bin);
        z_norm    = fliplr(z_norm);
    else
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = ensemble_profiles{nn}.re;
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = ensemble_profiles{nn}.re_CDP;
        end
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    re_bin(re_bin == 0)   = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0)   = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = z_norm >= bin_edges(bb) & z_norm <= bin_edges(bb+1);
        else
            idx_seg = z_norm >  bin_edges(bb) & z_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end

% Percentiles: [p10, p25, p50, p75, p90] per bin
pct = [10 50 90];
re_p_nd   = zeros(n_bins, length(pct));  re_p_d   = zeros(n_bins, length(pct));
lwc_p_nd  = zeros(n_bins, length(pct));  lwc_p_d  = zeros(n_bins, length(pct));
Nc_p_nd   = zeros(n_bins, length(pct));  Nc_p_d   = zeros(n_bins, length(pct));
vEff_p_nd = zeros(n_bins, length(pct));  vEff_p_d = zeros(n_bins, length(pct));

for bb = 1:n_bins
    re_p_nd(bb,:)   = prctile(vertSeg_nd{bb,1}, pct);
    lwc_p_nd(bb,:)  = prctile(vertSeg_nd{bb,2}, pct);
    Nc_p_nd(bb,:)   = prctile(vertSeg_nd{bb,3}, pct);
    vEff_p_nd(bb,:) = prctile(vertSeg_nd{bb,4}, pct);

    re_p_d(bb,:)    = prctile(vertSeg_d{bb,1}, pct);
    lwc_p_d(bb,:)   = prctile(vertSeg_d{bb,2}, pct);
    Nc_p_d(bb,:)    = prctile(vertSeg_d{bb,3}, pct);
    vEff_p_d(bb,:)  = prctile(vertSeg_d{bb,4}, pct);
end

disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])

% -------------------------------------------------------------------------
% ------ Plot VOCALS-REx panels (a)-(d) -----------------------------------
% -------------------------------------------------------------------------
figure1 = figure;

plt_clr_1   = mySavedColors(68,'fixed');
plt_clr_2   = mySavedColors(67,'fixed');

fnt_sz      = 18;
ttl_fnt     = 20;
alpha_inner = 0.30;
ln_w        = 1.5;
pnl_fnt     = 21;

y = [bin_center; flipud(bin_center)];

% --- (a) r_e ---
subplot(2,4,1)
hold on
% plot 10-90% data range as shaded region
fill([re_p_nd(:,1); flipud(re_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([re_p_d(:,1);  flipud(re_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 18])

annotation(figure1,'textbox',[0.130138 0.836 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (b) v_eff ---
subplot(2,4,2)
hold on
% plot 10-90% data range as shaded region
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

set(gca, 'XScale', 'log')
grid on; grid minor
xlim([0.005, 0.175])

annotation(figure1,'textbox',[0.3351315023418542 0.836 0.0373 0.0666],'String',{'(b)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (c) LWC ---
subplot(2,4,3)
hold on
% plot 10-90% data range as shaded region
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
annotation(figure1,'textbox',[0.540879493143727 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');
xlim([0, 0.725])

% --- (d) N_c ---
subplot(2,4,4)
hold on
% plot 10-90% data range as shaded region
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlim([0, 325])
annotation(figure1,'textbox',[0.746454846227638 0.836 0.0373 0.0666],'String','(d)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

legend('', 'Non-drizzling','', 'Drizzling', ...
    'Location','best','Interpreter','latex', ...
    'FontSize', fnt_sz-2, 'Color','white','TextColor','black', ...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])


% -------------------------------------------------------------------------
% ------ ORACLES: load, separate drizzle, compute percentiles -------------
% -------------------------------------------------------------------------
clearvars -except figure1 bin_center n_bins bin_edges fnt_sz ttl_fnt ...
    plt_clr_1 plt_clr_2 alpha_inner ln_w pnl_fnt

tauC_limit = 3;

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

end

load([foldername_data,...
    'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])

drizzle_LWP_threshold = 5;   % g/m^2
is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

vertSeg_nd = cell(n_bins, 4);
vertSeg_d  = cell(n_bins, 4);

num_removed = 0;
num_kept    = 0;

for nn = 1:N_profiles

    if max(ensemble_profiles{nn}.tau) < tauC_limit
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end

    % cloud profiles should range from cloud top to cloud bottom
    z_norm   = (ensemble_profiles{nn}.altitude - min(ensemble_profiles{nn}.altitude)) ./...
        (max(ensemble_profiles{nn}.altitude) - min(ensemble_profiles{nn}.altitude));


    significance_lvl = 0.1;
    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nd, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);
    v_eff_bin = 1 ./ (3 + gammaFit.alpha);

    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        re_bin    = fliplr(ensemble_profiles{nn}.re);
        lwc_bin   = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin    = fliplr(ensemble_profiles{nn}.total_Nc);
        v_eff_bin = flipud(v_eff_bin);
        z_norm    = fliplr(z_norm);
    else
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    re_bin(re_bin == 0)   = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0)   = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = z_norm >= bin_edges(bb) & z_norm <= bin_edges(bb+1);
        else
            idx_seg = z_norm >  bin_edges(bb) & z_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end

pct = [10 50 90];
re_p_nd   = zeros(n_bins, length(pct));  re_p_d   = zeros(n_bins, length(pct));
lwc_p_nd  = zeros(n_bins, length(pct));  lwc_p_d  = zeros(n_bins, length(pct));
Nc_p_nd   = zeros(n_bins, length(pct));  Nc_p_d   = zeros(n_bins, length(pct));
vEff_p_nd = zeros(n_bins, length(pct));  vEff_p_d = zeros(n_bins, length(pct));

for bb = 1:n_bins
    re_p_nd(bb,:)   = prctile(vertSeg_nd{bb,1}, pct);
    lwc_p_nd(bb,:)  = prctile(vertSeg_nd{bb,2}, pct);
    Nc_p_nd(bb,:)   = prctile(vertSeg_nd{bb,3}, pct);
    vEff_p_nd(bb,:) = prctile(vertSeg_nd{bb,4}, pct);

    re_p_d(bb,:)    = prctile(vertSeg_d{bb,1}, pct);
    lwc_p_d(bb,:)   = prctile(vertSeg_d{bb,2}, pct);
    Nc_p_d(bb,:)    = prctile(vertSeg_d{bb,3}, pct);
    vEff_p_d(bb,:)  = prctile(vertSeg_d{bb,4}, pct);
end

disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])

y = [bin_center; flipud(bin_center)];

% -------------------------------------------------------------------------
% ------ Plot ORACLES panels (e)-(h) --------------------------------------
% -------------------------------------------------------------------------

% --- (e) r_e ---
subplot(2,4,5)
hold on
% plot 10-90% data range as shaded region
fill([re_p_nd(:,1); flipud(re_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([re_p_d(:,1);  flipud(re_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(re_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Altitude', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 18])
annotation(figure1,'textbox',[0.132465746741153 0.367111111111112 0.0373 0.0666],'String',{'(e)'},...
    'Interpreter','latex','FontSize', pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (f) v_eff ---
subplot(2,4,6)
hold on
% plot 10-90% data range as shaded region
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(vEff_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

set(gca, 'XScale', 'log')
grid on; grid minor
xlabel('$\langle \nu_{e}(z) \rangle$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.005, 0.175])
annotation(figure1,'textbox',[0.343867473506009 0.367111111111112 0.0373 0.0666],'String',{'(f)'},...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (g) LWC ---
subplot(2,4,7)
hold on
% plot 10-90% data range as shaded region
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(lwc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.725])
annotation(figure1,'textbox',[0.542078635404322 0.361555555555556 0.0372999999999999 0.0666],'String','(g)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (h) N_c ---
subplot(2,4,8)
hold on
% plot 10-90% data range as shaded region
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,3))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_nd(:,2), bin_center, '-',  'Color', plt_clr_1, 'LineWidth', ln_w)
% plot 10-90% data range as shaded region
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,3))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
% plot the median value as a solid line
plot(Nc_p_d(:,2),  bin_center, '-',  'Color', plt_clr_2, 'LineWidth', ln_w)

grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 325])
annotation(figure1,'textbox',[0.745932819705434 0.367111111111112 0.0373 0.0666],'String','(h)',...
    'Interpreter','latex','FontSize',pnl_fnt,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% IEEE-column sizing
w = 7.16; h = 3;
figure1.Units    = 'inches';
figure1.Position = [1, 1, 2.75*w, 2.75*h];


% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle and non drizzle',...
    '- 4 panels with 10 and 90 percentiles.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle',...
    'and non drizzle - 8 panels with 10 and 90 percentiles - function of z.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------
