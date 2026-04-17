


%% Create median profiles of re, lwc and tau_c from ORACLES data

clear variables

% *** Optical Depth Filter ***
tauC_limit = 3;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

% keep track of the number of profiles kept and removed
num_removed = 0;
num_kept = 0;

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    % ** If tau_c is less than 3, skip this profile
    if max(tau_prof) < tauC_limit
        num_removed = num_removed +1;
        continue
    else
        num_kept = num_kept +1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd  = zeros(n_bins,1);  re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd = zeros(n_bins,1);  lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd  = zeros(n_bins,1);  Nc_iqr_nd  = zeros(n_bins,1);

re_med_d   = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d  = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d   = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_med_nd(bb)  = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb) = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)  = median(vertSeg_nd{bb,3}, 'omitnan');
    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_med_d(bb)   = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});
end



disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])


% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 22;
ttl_fnt = 26;


% Create axes
axes1 = axes('Parent',figure1,'Position',[0.185000000000001 0.11 0.21340579710145 0.815]);
hold(axes1,'on');

subplot(1,3,1)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 15])

% Create textbox
annotation(figure1,'textbox',[0.130138 0.806 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');









subplot(1,3,2)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.4])

% Create textbox
annotation(figure1,'textbox',[0.41393 0.806 0.0373 0.0666],'String','(b)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');







subplot(1,3,3)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 230])

% Create textbox
annotation(figure1,'textbox',[0.69073 0.806 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.871976270152657 0.882592592592596 0.122791322836865 0.111660079616977])


% Set the height and width as multiples of the IEEE paper requirements
w = 7.16; % inches
h = 3;    % inches

figure1.Units = 'inches';

figure1.Position = [1, 1, 2.5*w, 2.5*h]; % [left, bottom, width, height]




% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle',...
    '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------










%% Create median profiles of re, lwc and tau_c from VOCALS-REx data

clear variables

% *** Optical Depth Filter ***
tauC_limit = 3;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])
% load([foldername_data,...
%     'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25',...
%         '_drizzleLWP-threshold_5_10-Nov-2025.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

% keep track of the number of profiles kept and removed
num_removed = 0;
num_kept = 0;

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    % ** If tau_c is less than 3, skip this profile
    if max(tau_prof) < tauC_limit
        num_removed = num_removed +1;
        continue
    else
        num_kept = num_kept +1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re);
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re_CDP);
        end
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = ensemble_profiles{nn}.re;
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = ensemble_profiles{nn}.re_CDP;
        end
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd  = zeros(n_bins,1);  re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd = zeros(n_bins,1);  lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd  = zeros(n_bins,1);  Nc_iqr_nd  = zeros(n_bins,1);

re_med_d   = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d  = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d   = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_med_nd(bb)  = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb) = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)  = median(vertSeg_nd{bb,3}, 'omitnan');
    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_med_d(bb)   = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});
end




disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])



% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 22;
ttl_fnt = 26;


% Create axes
axes1 = axes('Parent',figure1,'Position',[0.185000000000001 0.11 0.21340579710145 0.815]);
hold(axes1,'on');

subplot(1,3,1)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)


% Create textbox
annotation(figure1,'textbox',[0.130138 0.806 0.0373 0.0666],'String',{'(d)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





subplot(1,3,2)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title('VOCALS-REx median profiles', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)
% Create textbox
annotation(figure1,'textbox',[0.41393 0.806 0.0373 0.0666],'String','(e)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');


subplot(1,3,3)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% Create textbox
annotation(figure1,'textbox',[0.69073 0.806 0.0373 0.0666],'String','(f)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.871976270152657 0.882592592592596 0.122791322836865 0.111660079616977])




% Set the height and width as multiples of the IEEE paper requirements
w = 7.16; % inches
h = 3;    % inches

figure1.Units = 'inches';

figure1.Position = [1, 1, 2.5*w, 2.5*h]; % [left, bottom, width, height]



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle',...
    '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------

















%% Combine VOCALS-REx and ORACLES in a subplot to show median profiles and the IQR 


% -------------------------------------------------------------------------
% ------ Create median profiles of re, lwc and tau_c from VOCALS-REx data ---
% -------------------------------------------------------------------------
clear variables

% *** Optical Depth Filter ***
tauC_limit = 3;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])
% load([foldername_data,...
%     'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25',...
%         '_drizzleLWP-threshold_5_10-Nov-2025.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

% keep track of the number of profiles kept and removed
num_removed = 0;
num_kept = 0;

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    % ** If tau_c is less than 3, skip this profile
    if max(tau_prof) < tauC_limit
        num_removed = num_removed +1;
        continue
    else
        num_kept = num_kept +1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re);
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = fliplr(ensemble_profiles{nn}.re_CDP);
        end
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        if isfield(ensemble_profiles{nn}, 're')==true
            re_bin  = ensemble_profiles{nn}.re;
        elseif isfield(ensemble_profiles{nn}, 're_CDP')==true
            re_bin  = ensemble_profiles{nn}.re_CDP;
        end
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd  = zeros(n_bins,1);  re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd = zeros(n_bins,1);  lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd  = zeros(n_bins,1);  Nc_iqr_nd  = zeros(n_bins,1);

re_med_d   = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d  = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d   = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_med_nd(bb)  = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb) = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)  = median(vertSeg_nd{bb,3}, 'omitnan');
    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_med_d(bb)   = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});
end




disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])



% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 18;
ttl_fnt = 20;


% Create axes
axes1 = axes('Parent',figure1,'Position',[0.185000000000001 0.11 0.21340579710145 0.815]);
hold(axes1,'on');

subplot(2,3,1)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 15])

% Create textbox
annotation(figure1,'textbox',[0.130138 0.836 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





subplot(2,3,2)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['VOCALS-REx median profiles - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
% Create textbox
annotation(figure1,'textbox',[0.41393 0.836 0.0373 0.0666],'String','(b)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 0.675])

subplot(2,3,3)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% Create textbox
annotation(figure1,'textbox',[0.69073 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 240])

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])




% Set the height and width as multiples of the IEEE paper requirements
w = 7.16; % inches
h = 3;    % inches

figure1.Units = 'inches';

figure1.Position = [1, 1, 2.5*w, 2.5*h]; % [left, bottom, width, height]










% -------------------------------------------------------------------------
% ------ Create median profiles of re, lwc and tau_c from ORACLES data ---
% -------------------------------------------------------------------------

clearvars -except figure1

% *** Optical Depth Filter ***
tauC_limit = 3;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')

    foldername_data = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Hyperspectral_Cloud_Retrievals/ORACLES/oracles_data/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

% keep track of the number of profiles kept and removed
num_removed = 0;
num_kept = 0;

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    % ** If tau_c is less than 3, skip this profile
    if max(tau_prof) < tauC_limit
        num_removed = num_removed +1;
        continue
    else
        num_kept = num_kept +1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
    else
        % Descending: data already starts at cloud top
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0) = 0.01;
    lwc_bin(lwc_bin == 0) = 0.001;
    Nc_bin(Nc_bin == 0) = 0.01;

    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd  = zeros(n_bins,1);  re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd = zeros(n_bins,1);  lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd  = zeros(n_bins,1);  Nc_iqr_nd  = zeros(n_bins,1);

re_med_d   = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d  = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d   = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    re_med_nd(bb)  = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb) = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)  = median(vertSeg_nd{bb,3}, 'omitnan');
    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    re_med_d(bb)   = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)  = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)   = median(vertSeg_d{bb,3}, 'omitnan');
    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});
end



disp(['# profiles kept   : ', num2str(num_kept)])
disp(['# profiles removed: ', num2str(num_removed)])


% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz = 18;
ttl_fnt = 20;


% Create axes
axes1 = axes('Parent',figure1,'Position',[0.185000000000001 0.11 0.21340579710145 0.815]);
hold(axes1,'on');

subplot(2,3,4)
hold on
x = [re_med_nd - re_iqr_nd/2; flipud(re_med_nd + re_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_iqr_d/2; flipud(re_med_d + re_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(z) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 15])

% Create textbox
annotation(figure1,'textbox',[0.132465746741155 0.367111111111112 0.0373 0.0666],'String',{'(d)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');









subplot(2,3,5)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_iqr_nd/2; flipud(lwc_med_nd + lwc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_iqr_d/2; flipud(lwc_med_d + lwc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(z) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.675])

% Create textbox
annotation(figure1,'textbox',[0.417226764742396 0.361555555555556 0.0372999999999999 0.0666],'String','(e)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');







subplot(2,3,6)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_iqr_nd/2; flipud(Nc_med_nd + Nc_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_iqr_d/2; flipud(Nc_med_d + Nc_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(z) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 240])

% Create textbox
annotation(figure1,'textbox',[0.688797217877096 0.367111111111112 0.0373 0.0666],'String','(f)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
end
saveas(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle and non drizzle.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle and non drizzle',...
    '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------














%% Create imagesc figures of the droplet distribution from VOCALS-REx with overlaid median profile

clear variables


three_panel_flag = false;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/NCAR_C130/SPS_1/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_25_04-Dec-2025.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DC >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc, full n(r)}
vertSeg_d  = cell(n_bins, 4);

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        if ensemble_profiles{nn}.flag_2DC_data_is_conforming == false

            re_bin  = fliplr(ensemble_profiles{nn}.re_CDP);

        else

            re_bin  = fliplr(ensemble_profiles{nn}.re);

        end

        lwc_bin     = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin      = fliplr(ensemble_profiles{nn}.total_Nc);
        % we also need to account for the Nc distribution being less
        % than the typical size
        Nc_dist_bin = zeros(91, size(ensemble_profiles{nn}.Nc, 2));
        Nc_dist_bin(1:size(ensemble_profiles{nn}.Nc, 1), :) = fliplr(ensemble_profiles{nn}.Nc);

    else
        % Descending: data already starts at cloud top
        if ensemble_profiles{nn}.flag_2DC_data_is_conforming == false

            re_bin  = ensemble_profiles{nn}.re_CDP;


        else

            re_bin  = ensemble_profiles{nn}.re;


        end

        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
        % we also need to account for the Nc distribution being less
        % than the typical size
        Nc_dist_bin = zeros(91, size(ensemble_profiles{nn}.Nc, 2));
        Nc_dist_bin(1:size(ensemble_profiles{nn}.Nc, 1), :) = ensemble_profiles{nn}.Nc;

    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0)           = 0.01;
    lwc_bin(lwc_bin == 0)         = 0.001;
    Nc_bin(Nc_bin == 0)           = 0.01;
    Nc_dist_bin(Nc_dist_bin == 0) = 0.0001;


    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  Nc_dist_bin(:, idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, Nc_dist_bin(:, idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd       = zeros(n_bins,1);   re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd      = zeros(n_bins,1);   lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd       = zeros(n_bins,1);   Nc_iqr_nd  = zeros(n_bins,1);
Nc_dist_med_nd  = zeros(91, n_bins);
Nc_dist_mean_nd = zeros(91, n_bins);

re_med_d       = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d      = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d       = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);
Nc_dist_med_d  = zeros(91, n_bins);
Nc_dist_mean_d = zeros(91, n_bins);

Nc_dist_med_all  = zeros(91, n_bins);
Nc_dist_mean_all = zeros(91, n_bins);
re_med_all       = zeros(n_bins,1);  re_iqr_all   = zeros(n_bins,1); re_std_all = zeros(n_bins, 1);


for bb = 1:n_bins
    
    % --- profiles without drizzle ---
    re_med_nd(bb)          = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)         = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)          = median(vertSeg_nd{bb,3}, 'omitnan');
    Nc_dist_med_nd(:,bb)   = median(vertSeg_nd{bb,4}, 2, 'omitnan');
    Nc_dist_mean_nd(:,bb)  = mean(vertSeg_nd{bb,4}, 2, 'omitnan');

    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    % --- profiles with drizzle ---
    re_med_d(bb)           = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)          = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)           = median(vertSeg_d{bb,3}, 'omitnan');
    Nc_dist_med_nd(:,bb)   = median(vertSeg_nd{bb,4}, 2, 'omitnan');
    Nc_dist_mean_nd(:,bb)  = mean(vertSeg_nd{bb,4}, 2, 'omitnan');

    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});

    % --- All Profiles ----
    Nc_dist_med_all(:, bb)  = median([vertSeg_nd{bb,4}, vertSeg_d{bb,4}], 2, 'omitnan');
    Nc_dist_mean_all(:, bb) = mean([vertSeg_nd{bb,4}, vertSeg_d{bb,4}], 2, 'omitnan');

    re_med_all(bb)          = median([vertSeg_nd{bb,1}, vertSeg_d{bb,1}], 'omitnan');
    re_iqr_all(bb)          = iqr([vertSeg_nd{bb,1}, vertSeg_d{bb,1}]);
    re_std_all(bb)          = std([vertSeg_nd{bb,1}, vertSeg_d{bb,1}]); 

end






% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz  = 22;
ttl_fnt = 26;
ax_fnt  = 22;
cb_fnt  = 20;

% -------------------------------------------------------------------------
% Extract data
% -------------------------------------------------------------------------
norm_optical_depth   = bin_center;                                   % 1x28, unitless
r_bins_all = double(ensemble_profiles{1}.drop_radius_bin_center);    % 1x91, microns

% plot bins less than or equal to 50 microns
idx_lessThan50 = r_bins_all <= 50;


% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_mean_nd);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax1 = subplot(1, 3, 1);
else
    ax1 = axes;

end


imagesc(ax1, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax1, 'YDir', 'reverse');     % altitude increases upward
colormap(ax1, 'turbo');

cb1 = colorbar(ax1);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax1, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax1, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax1, ['Mean $N_{c} $ distribution - VOCALS-REx', newline, 'No Drizzle'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])





figure2 = figure;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_med_nd);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax2 = subplot(1, 3, 1);
else
    ax2 = axes;

end


imagesc(ax2, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax2, 'YDir', 'reverse');     % altitude increases upward
colormap(ax2, 'turbo');

cb1 = colorbar(ax2);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax2, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax2, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax2, ['Median $N_{c} $ distribution - VOCALS-REx', newline, 'No Drizzle'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])





% --- Plot the median of all profiles ---
figure3 = figure;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_med_all);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax3 = subplot(1, 3, 1);
else
    ax3 = axes;

end


imagesc(ax3, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax3, 'YDir', 'reverse');     % altitude increases upward
colormap(ax3, 'turbo');

cb1 = colorbar(ax3);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax3, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax3, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax3, ['Median $N_{c} $ distribution - VOCALS-REx', newline, 'All Profiles'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])






% --- Plot the mean of all profiles ---
figure4 = figure;
ln_sz = 2.5;
lgnd_fnt = 20;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_mean_all);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax4 = subplot(1, 3, 1);
else
    ax4 = axes;

end


imagesc(ax4, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax4, 'YDir', 'reverse');     % altitude increases upward
colormap(ax4, 'turbo');


% overlay median r_eff profile and IQR lines
hold(ax4, 'on')
plot(ax4, re_med_all, norm_optical_depth, 'k-',  'LineWidth', ln_sz,   'DisplayName', 'Median $r_{e}$')
% plot(ax4, re_med_all - re_iqr_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'DisplayName', '$\mu \pm \sigma$')
% plot(ax4, re_med_all + re_iqr_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'HandleVisibility', 'off')
plot(ax4, re_med_all - re_std_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'DisplayName', '$\mu \pm \sigma$')
plot(ax4, re_med_all + re_std_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'HandleVisibility', 'off')
hold(ax4, 'off')

legend(ax4, 'Interpreter', 'latex', 'Location', 'best', ...
    'Color', 'w', 'TextColor', 'k', 'FontSize', lgnd_fnt)


cb1 = colorbar(ax4);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax4, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax4, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax4, ['Mean $N_{c} $ distribution - VOCALS-REx', newline, 'All Profiles'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])




% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
% end
% saveas(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle.fig']);
%
%
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle',...
%     '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------







%% Create imagesc figures of the droplet distribution from ORACLES with overlaid median profile

clear variables


three_panel_flag = false;

% *** Optical Depth Filter ***
idx_tauC = 3;

% -------------------------------------------------------------------------
%  File locations
%  -------------------------------------------------------------------------

which_computer = whatComputer;

if strcmp(which_computer, 'anbu8374')



elseif strcmp(which_computer, 'andrewbuggee')

    foldername_data = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/ORACLES/ORACLES_data/'];

    foldername_save = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/', ...
        'Presentations_and_Papers/paper_3/'];

end


% -------------------------------------------------------------------------
%  Load ensemble profiles mat file
%  -------------------------------------------------------------------------

load([foldername_data,...
    'ensemble_profiles_with_precip_from_33_files_LWC-threshold_0.05_Nc-threshold_10_no_rEff_greaterThan50_microns_16-Mar-2026.mat'])

N_profiles = length(ensemble_profiles);
disp(['Loaded ', num2str(N_profiles), ' vertical profiles.'])




% -------------------------------------------------------------------------
%  1. Separate profiles: drizzling vs non-drizzling
%     Criterion: rain/drizzle liquid water path >= threshold
%  -------------------------------------------------------------------------

drizzle_LWP_threshold = 5;   % g/m²

is_drizzle = false(1, N_profiles);
for nn = 1:N_profiles
    is_drizzle(nn) = ensemble_profiles{nn}.lwp_2DS_HVPS >= drizzle_LWP_threshold;
end

idx_drizzle    = find(is_drizzle);
idx_no_drizzle = find(~is_drizzle);

disp(['Drizzling profiles:     ', num2str(length(idx_drizzle))])
disp(['Non-drizzling profiles: ', num2str(length(idx_no_drizzle))])





% -------------------------------------------------------------------------
%      Ensemble MEDIAN profiles of r_e, LWC, N_c vs normalized optical depth
%         for non-precipitating (no-drizzle) clouds
%
%  Analogous to first_paper_figures.m (VOCALS-REx version). Take 1
%  -------------------------------------------------------------------------

n_bins    = 30;
bin_edges = 0:1/n_bins:1;
bin_center = ((bin_edges(1:end-1) + bin_edges(2:end)) / 2)';

% Accumulate by vertical bin (normalized optical depth, 0 = cloud top, 1 = base)
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc, full n(r)}
vertSeg_d  = cell(n_bins, 4);

for nn = 1:N_profiles

    % Normalize optical depth to [0,1] (0 = cloud top, 1 = cloud base)
    tau_prof = ensemble_profiles{nn}.tau;
    tau_norm = tau_prof ./ max(tau_prof);

    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0

        % Ascending: flip so cloud top is first
        re_bin      = fliplr(ensemble_profiles{nn}.re);
        lwc_bin     = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin      = fliplr(ensemble_profiles{nn}.total_Nc);
        % we also need to account for the Nc distribution being less
        % than the typical size
        Nc_dist_bin = zeros(172, size(ensemble_profiles{nn}.Nd, 2));

        if size(ensemble_profiles{nn}.Nd,1) == 172

            Nc_dist_bin(1:size(ensemble_profiles{nn}.Nd, 1), :) = fliplr(ensemble_profiles{nn}.Nd);

        elseif size(ensemble_profiles{nn}.Nd,1) > 172

            Nc_dist_bin(1:172, :) = fliplr(ensemble_profiles{nn}.Nd(1:172, :));
        end

    else

        % Descending: data already starts at cloud top
        re_bin  = ensemble_profiles{nn}.re;
        lwc_bin = ensemble_profiles{nn}.lwc;
        Nc_bin  = ensemble_profiles{nn}.total_Nc;
        % we also need to account for the Nc distribution being less
        % than the typical size
        Nc_dist_bin = zeros(172, size(ensemble_profiles{nn}.Nd, 2));

        if size(ensemble_profiles{nn}.Nd,1) == 172
            Nc_dist_bin(1:size(ensemble_profiles{nn}.Nd, 1), :) = ensemble_profiles{nn}.Nd;

        elseif size(ensemble_profiles{nn}.Nd,1) > 172
            Nc_dist_bin(1:172, :) = ensemble_profiles{nn}.Nd(1:172, :);
        end

    end

    % get rid of zeros
    % Remove zeros from the accumulated data
    re_bin(re_bin == 0)           = 0.01;
    lwc_bin(lwc_bin == 0)         = 0.001;
    Nc_bin(Nc_bin == 0)           = 0.01;
    Nc_dist_bin(Nc_dist_bin == 0) = 0.0001;


    for bb = 1:n_bins
        if bb == 1
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        end

        if is_drizzle(nn)
            vertSeg_d{bb,1} = [vertSeg_d{bb,1},  re_bin(idx_seg)];
            vertSeg_d{bb,2} = [vertSeg_d{bb,2},  lwc_bin(idx_seg)];
            vertSeg_d{bb,3} = [vertSeg_d{bb,3},  Nc_bin(idx_seg)];
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  Nc_dist_bin(:, idx_seg)];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, Nc_dist_bin(:, idx_seg)];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd       = zeros(n_bins,1);   re_iqr_nd  = zeros(n_bins,1);
lwc_med_nd      = zeros(n_bins,1);   lwc_iqr_nd = zeros(n_bins,1);
Nc_med_nd       = zeros(n_bins,1);   Nc_iqr_nd  = zeros(n_bins,1);
Nc_dist_med_nd  = zeros(172, n_bins);
Nc_dist_mean_nd = zeros(172, n_bins);

re_med_d       = zeros(n_bins,1);  re_iqr_d   = zeros(n_bins,1);
lwc_med_d      = zeros(n_bins,1);  lwc_iqr_d  = zeros(n_bins,1);
Nc_med_d       = zeros(n_bins,1);  Nc_iqr_d   = zeros(n_bins,1);
Nc_dist_med_d  = zeros(172, n_bins);
Nc_dist_mean_d = zeros(172, n_bins);

Nc_dist_med_all  = zeros(172, n_bins);
Nc_dist_mean_all = zeros(172, n_bins);
re_med_all       = zeros(n_bins,1);  re_iqr_all   = zeros(n_bins,1); re_std_all = zeros(n_bins, 1);


for bb = 1:n_bins
    
    % --- profiles without drizzle ---
    re_med_nd(bb)          = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)         = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)          = median(vertSeg_nd{bb,3}, 'omitnan');
    Nc_dist_med_nd(:,bb)   = median(vertSeg_nd{bb,4}, 2, 'omitnan');
    Nc_dist_mean_nd(:,bb)  = mean(vertSeg_nd{bb,4}, 2, 'omitnan');

    re_iqr_nd(bb)  = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb) = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)  = iqr(vertSeg_nd{bb,3});

    % --- profiles with drizzle ---
    re_med_d(bb)           = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)          = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)           = median(vertSeg_d{bb,3}, 'omitnan');
    Nc_dist_med_nd(:,bb)   = median(vertSeg_nd{bb,4}, 2, 'omitnan');
    Nc_dist_mean_nd(:,bb)  = mean(vertSeg_nd{bb,4}, 2, 'omitnan');

    re_iqr_d(bb)   = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)  = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)   = iqr(vertSeg_d{bb,3});

    % --- All Profiles ----
    Nc_dist_med_all(:, bb)  = median([vertSeg_nd{bb,4}, vertSeg_d{bb,4}], 2, 'omitnan');
    Nc_dist_mean_all(:, bb) = mean([vertSeg_nd{bb,4}, vertSeg_d{bb,4}], 2, 'omitnan');

    re_med_all(bb)          = median([vertSeg_nd{bb,1}, vertSeg_d{bb,1}], 'omitnan');
    re_iqr_all(bb)          = iqr([vertSeg_nd{bb,1}, vertSeg_d{bb,1}]);
    re_std_all(bb)          = std([vertSeg_nd{bb,1}, vertSeg_d{bb,1}]); 

end






% --- Overlay drizzling vs non-drizzling on same axes ---
% Create figure
figure1 = figure;

plt_clr_1 = mySavedColors(64,'fixed');
plt_clr_2 = mySavedColors(62,'fixed');

fnt_sz  = 22;
ttl_fnt = 26;
ax_fnt  = 22;
cb_fnt  = 20;

% -------------------------------------------------------------------------
% Extract data
% -------------------------------------------------------------------------
norm_optical_depth   = bin_center;                                   % 1x28, unitless
r_bins_all = double(ensemble_profiles{1}.drop_radius_bin_center);    % 1x91, microns

% plot bins less than or equal to 50 microns
idx_lessThan50 = r_bins_all <= 50;


% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_mean_nd);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax1 = subplot(1, 3, 1);
else
    ax1 = axes;

end


imagesc(ax1, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax1, 'YDir', 'reverse');     % altitude increases upward
colormap(ax1, 'turbo');

cb1 = colorbar(ax1);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax1, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax1, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax1, ['Mean $N_{c} $ distribution - ORACLES', newline, 'No Drizzle'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])





figure2 = figure;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_med_nd);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax2 = subplot(1, 3, 1);
else
    ax2 = axes;

end


imagesc(ax2, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax2, 'YDir', 'reverse');     % altitude increases upward
colormap(ax2, 'turbo');

cb1 = colorbar(ax2);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax2, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax2, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax2, ['Median $N_{c} $ distribution - ORACLES', newline, 'No Drizzle'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])





% --- Plot the median of all profiles ---
figure3 = figure;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_med_all);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax3 = subplot(1, 3, 1);
else
    ax3 = axes;

end


imagesc(ax3, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax3, 'YDir', 'reverse');     % altitude increases upward
colormap(ax3, 'turbo');

cb1 = colorbar(ax3);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax3, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax3, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax3, ['Median $N_{c} $ distribution - ORACLES', newline, 'All Profiles'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])






% --- Plot the mean of all profiles ---
figure4 = figure;
ln_sz = 2.5;
lgnd_fnt = 20;

% -------------------------------------------------------------------------
% Prepare log10 of Nc for display (zeros/NaNs -> NaN so they show as background)
% -------------------------------------------------------------------------
Nc_log = log10(Nc_dist_mean_all);
Nc_log(~isfinite(Nc_log)) = NaN;

% imagesc expects the matrix as [n_y x n_x] = [n_tau x n_bins]
Nc_img = Nc_log';    % n_tau x n_bins


if three_panel_flag
    ax4 = subplot(1, 3, 1);
else
    ax4 = axes;

end


imagesc(ax4, r_bins_all(idx_lessThan50), norm_optical_depth, Nc_img(:, idx_lessThan50));
set(ax4, 'YDir', 'reverse');     % altitude increases upward
colormap(ax4, 'turbo');


% overlay median r_eff profile and IQR lines
hold(ax4, 'on')
plot(ax4, re_med_all, norm_optical_depth, 'k-',  'LineWidth', ln_sz,   'DisplayName', 'Median $r_{e}$')
% plot(ax4, re_med_all - re_iqr_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'DisplayName', '$\mu \pm \sigma$')
% plot(ax4, re_med_all + re_iqr_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'HandleVisibility', 'off')
plot(ax4, re_med_all - re_std_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'DisplayName', '$\mu \pm \sigma$')
plot(ax4, re_med_all + re_std_all,   norm_optical_depth, 'k--', 'LineWidth', ln_sz-0.5, 'HandleVisibility', 'off')
hold(ax4, 'off')

legend(ax4, 'Interpreter', 'latex', 'Location', 'best', ...
    'Color', 'w', 'TextColor', 'k', 'FontSize', lgnd_fnt)


cb1 = colorbar(ax4);
cb1.Label.String     = '$\log_{10}(N_c \ [\mathrm{cm}^{-3}])$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize   = cb_fnt;

xlabel(ax4, 'Droplet radius ($\mu$m)', 'Interpreter', 'latex', 'FontSize', ax_fnt)
ylabel(ax4, 'Normalized Optical Depth',            'Interpreter', 'latex', 'FontSize', ax_fnt)
title(ax4, ['Mean $N_{c} $ distribution - ORACLES', newline, 'All Profiles'], ...
    'Interpreter', 'latex', 'FontSize', ttl_fnt)

set(gcf, 'Position', [0 0 1200 600])




% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
% end
% saveas(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle.fig']);
%
%
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(figure1,[folderpath_figs,'VOCALS-REx in-situ median profiles separated by drizzle and non drizzle',...
%     '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------

