


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
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
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

figure1.Position = [1, 1, 2.75*w, 2.75*h]; % [left, bottom, width, height]










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
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
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
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
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


















%% Combine VOCALS-REx and ORACLES in a subplot to show median profiles and the IQR 
% *** 4 panels including effective variance ***
% *** Plot the Inter-quartile range (25% - 75%) as the uncertainty


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
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 4);

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



    % ------------------------------------------------------------------
    % Fit the droplet size distribution to get the gamma shape parameter
    % This is the ORACLES equivalent of the pre-stored gammaFit.alpha from
    % VOCALS-REx profiles.  We compute it here because find_verticalProfiles_ORACLES
    % does not store the fit result on the profile struct.

    significance_lvl = 0.1;

    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nc, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);

    v_eff_bin = 1 ./ (3 + gammaFit.alpha);            % per-altitude effective variance
    % -------------------------------------------------------------------------




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
        v_eff_bin = flipud(v_eff_bin);
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
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd    = zeros(n_bins,1);  re_iqr_nd    = zeros(n_bins,1);
lwc_med_nd   = zeros(n_bins,1);  lwc_iqr_nd   = zeros(n_bins,1);
Nc_med_nd    = zeros(n_bins,1);  Nc_iqr_nd    = zeros(n_bins,1);
vEff_med_nd  = zeros(n_bins,1);  vEff_iqr_nd  = zeros(n_bins,1);

re_med_d     = zeros(n_bins,1);  re_iqr_d     = zeros(n_bins,1);
lwc_med_d    = zeros(n_bins,1);  lwc_iqr_d    = zeros(n_bins,1);
Nc_med_d     = zeros(n_bins,1);  Nc_iqr_d     = zeros(n_bins,1);
vEff_med_d   = zeros(n_bins,1);  vEff_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins

     % --- without drizzle ---
    re_med_nd(bb)    = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)   = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)    = median(vertSeg_nd{bb,3}, 'omitnan');
    vEff_med_nd(bb)  = median(vertSeg_nd{bb,4}, 'omitnan');

    re_iqr_nd(bb)    = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb)   = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)    = iqr(vertSeg_nd{bb,3});
    vEff_iqr_nd(bb)  = iqr(vertSeg_nd{bb,4});


    % --- with drizzle ---
    re_med_d(bb)     = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)    = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)     = median(vertSeg_d{bb,3}, 'omitnan');
    vEff_med_d(bb)   = median(vertSeg_d{bb,4}, 'omitnan');
    
    re_iqr_d(bb)    = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)   = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)    = iqr(vertSeg_d{bb,3});
    vEff_iqr_d(bb)  = iqr(vertSeg_d{bb,4});

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

subplot(2,4,1)
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






subplot(2,4,2)
hold on
x = [vEff_med_nd - vEff_iqr_nd/2; flipud(vEff_med_nd + vEff_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [vEff_med_d - vEff_iqr_d/2; flipud(vEff_med_d + vEff_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.005, 0.175])

% Create textbox
annotation(figure1,'textbox',[0.451315023418542 0.836 0.0373 0.0666],'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');








subplot(2,4,3)

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
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
% Create textbox
annotation(figure1,'textbox',[0.540879493143727 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 0.675])






subplot(2,4,4)

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
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% Create textbox
annotation(figure1,'textbox',[0.746454846227638 0.836 0.0373 0.0666],'String','(d)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 240])

legend({'','Non-drizzling',...
    '','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])
% legend({'',['Non-drizzling - N = ', num2str(length(idx_no_drizzle))],...
%     '',['Drizzling - N = ', num2str(length(idx_drizzle))]}, 'Location', 'best',...
%     'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
%     'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])




% Set the height and width as multiples of the IEEE paper requirements
w = 7.16; % inches
h = 3;    % inches

figure1.Units = 'inches';

figure1.Position = [1, 1, 2.75*w, 2.75*h]; % [left, bottom, width, height]










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
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc, var_eff}
vertSeg_d  = cell(n_bins, 4);

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


    % ------------------------------------------------------------------
    % Fit the droplet size distribution to get the gamma shape parameter
    % This is the ORACLES equivalent of the pre-stored gammaFit.alpha from
    % VOCALS-REx profiles.  We compute it here because find_verticalProfiles_ORACLES
    % does not store the fit result on the profile struct.

    significance_lvl = 0.1;

    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nd, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);

    v_eff_bin = 1 ./ (3 + gammaFit.alpha);            % per-altitude effective variance
    % -------------------------------------------------------------------------


    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
        v_eff_bin = flipud(v_eff_bin);
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
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end

    end

end




% Compute median and IQR for each bin
re_med_nd    = zeros(n_bins,1);  re_iqr_nd    = zeros(n_bins,1);
lwc_med_nd   = zeros(n_bins,1);  lwc_iqr_nd   = zeros(n_bins,1);
Nc_med_nd    = zeros(n_bins,1);  Nc_iqr_nd    = zeros(n_bins,1);
vEff_med_nd  = zeros(n_bins,1);  vEff_iqr_nd  = zeros(n_bins,1);

re_med_d     = zeros(n_bins,1);  re_iqr_d     = zeros(n_bins,1);
lwc_med_d    = zeros(n_bins,1);  lwc_iqr_d    = zeros(n_bins,1);
Nc_med_d     = zeros(n_bins,1);  Nc_iqr_d     = zeros(n_bins,1);
vEff_med_d   = zeros(n_bins,1);  vEff_iqr_d   = zeros(n_bins,1);

for bb = 1:n_bins
    % --- without drizzle ---
    re_med_nd(bb)    = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)   = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)    = median(vertSeg_nd{bb,3}, 'omitnan');
    vEff_med_nd(bb)  = median(vertSeg_nd{bb,4}, 'omitnan');

    re_iqr_nd(bb)    = iqr(vertSeg_nd{bb,1});
    lwc_iqr_nd(bb)   = iqr(vertSeg_nd{bb,2});
    Nc_iqr_nd(bb)    = iqr(vertSeg_nd{bb,3});
    vEff_iqr_nd(bb)  = iqr(vertSeg_nd{bb,4});


    % --- with drizzle ---
    re_med_d(bb)     = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)    = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)     = median(vertSeg_d{bb,3}, 'omitnan');
    vEff_med_d(bb)   = median(vertSeg_d{bb,4}, 'omitnan');
    
    re_iqr_d(bb)    = iqr(vertSeg_d{bb,1});
    lwc_iqr_d(bb)   = iqr(vertSeg_d{bb,2});
    Nc_iqr_d(bb)    = iqr(vertSeg_d{bb,3});
    vEff_iqr_d(bb)  = iqr(vertSeg_d{bb,4});

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

subplot(2,4,5)
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
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 15])

% Create textbox
annotation(figure1,'textbox',[0.132465746741153 0.367111111111112 0.0373 0.0666],'String',{'(e)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





subplot(2,4,6)
hold on
x = [vEff_med_nd - vEff_iqr_nd/2; flipud(vEff_med_nd + vEff_iqr_nd/2)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [vEff_med_d - vEff_iqr_d/2; flipud(vEff_med_d + vEff_iqr_d/2)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
xlabel('$\langle \nu_{e}(\tau) \rangle$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.005, 0.175])

% Create textbox
annotation(figure1,'textbox',[0.343867473506009 0.367111111111112 0.0373 0.0666],'String',{'(f)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');









subplot(2,4,7)

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
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.675])

% Create textbox
annotation(figure1,'textbox',[0.542078635404322 0.361555555555556 0.0372999999999999 0.0666],'String','(g)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');







subplot(2,4,8)

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
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 240])

% Create textbox
annotation(figure1,'textbox',[0.745932819705434 0.367111111111112 0.0373 0.0666],'String','(h)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% legend({'',['Non-drizzling - N = ', num2str(length(idx_no_drizzle))],...
%     '',['Drizzling - N = ', num2str(length(idx_drizzle))]}, 'Location', 'best',...
%     'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
%     'Position',[0.871976270152657 0.367111111111112 0.122791322836865 0.111660079616977])





% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_3/saved_figures/';
% end
% saveas(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle and non drizzle',...
%     '- 4 panels with IQR.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle',...
%     'and non drizzle - 4 panels with IQR.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------















%% Combine VOCALS-REx and ORACLES in a subplot to show median profiles and the IQR 
% *** 4 panels including effective variance ***
% *** With asymmetric uncert ***

% I don't know if the Interquartile range is the best way to display the
% spread in the values of each variable. These variables tend to follow log
% normal distributions, and the IQR doesn't show the asymmetric nature of
% the distribution. Try plotting the average deviation above and below the
% median value at each level instead. 


% ------------------------------------------------------------------------
% ------ Create median profiles of re, vEff, lwc, and tau_c from VOCALS-REx data ---
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
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 4);

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



    % ------------------------------------------------------------------
    % Fit the droplet size distribution to get the gamma shape parameter
    % This is the ORACLES equivalent of the pre-stored gammaFit.alpha from
    % VOCALS-REx profiles.  We compute it here because find_verticalProfiles_ORACLES
    % does not store the fit result on the profile struct.

    significance_lvl = 0.1;

    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nc, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);

    v_eff_bin = 1 ./ (3 + gammaFit.alpha);            % per-altitude effective variance
    % -------------------------------------------------------------------------




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
        v_eff_bin = flipud(v_eff_bin);
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
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end
    end

end




% Compute median and IQR for each bin
re_med_nd    = zeros(n_bins,1);  re_avgDevAbove_nd    = zeros(n_bins,1); re_avgDevBelow_nd    = zeros(n_bins,1);
lwc_med_nd   = zeros(n_bins,1);  lwc_avgDevAbove_nd   = zeros(n_bins,1); lwc_avgDevBelow_nd   = zeros(n_bins,1);
Nc_med_nd    = zeros(n_bins,1);  Nc_avgDevAbove_nd    = zeros(n_bins,1); Nc_avgDevBelow_nd    = zeros(n_bins,1);
vEff_med_nd  = zeros(n_bins,1);  vEff_avgDevAbove_nd  = zeros(n_bins,1); vEff_avgDevBelow_nd  = zeros(n_bins,1);

re_med_d     = zeros(n_bins,1);  re_avgDevAbove_d    = zeros(n_bins,1); re_avgDevBelow_d    = zeros(n_bins,1);
lwc_med_d    = zeros(n_bins,1);  lwc_avgDevAbove_d   = zeros(n_bins,1); lwc_avgDevBelow_d   = zeros(n_bins,1);
Nc_med_d     = zeros(n_bins,1);  Nc_avgDevAbove_d    = zeros(n_bins,1); Nc_avgDevBelow_d    = zeros(n_bins,1);
vEff_med_d   = zeros(n_bins,1);  vEff_avgDevAbove_d  = zeros(n_bins,1); vEff_avgDevBelow_d  = zeros(n_bins,1);

for bb = 1:n_bins

     % --- without drizzle ---
    re_med_nd(bb)    = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)   = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)    = median(vertSeg_nd{bb,3}, 'omitnan');
    vEff_med_nd(bb)  = median(vertSeg_nd{bb,4}, 'omitnan');

    re = vertSeg_nd{bb,1};
    re_avgDevAbove_nd(bb)    = mean( re(re > re_med_nd(bb)) - re_med_nd(bb));
    re_avgDevBelow_nd(bb)    = mean( re_med_nd(bb) - re(re < re_med_nd(bb)));

    lwc = vertSeg_nd{bb,2};
    lwc_avgDevAbove_nd(bb)    = mean( lwc(lwc > lwc_med_nd(bb)) - lwc_med_nd(bb));
    lwc_avgDevBelow_nd(bb)    = mean( lwc_med_nd(bb) - lwc(lwc < lwc_med_nd(bb)));

    Nc = vertSeg_nd{bb,3};
    Nc_avgDevAbove_nd(bb)    = mean( Nc(Nc > Nc_med_nd(bb)) - Nc_med_nd(bb));
    Nc_avgDevBelow_nd(bb)    = mean( Nc_med_nd(bb) - Nc(Nc < Nc_med_nd(bb)));

    vEff = vertSeg_nd{bb,4};
    vEff_avgDevAbove_nd(bb)    = mean( vEff(vEff > vEff_med_nd(bb)) - vEff_med_nd(bb));
    vEff_avgDevBelow_nd(bb)    = mean( vEff_med_nd(bb) - vEff(vEff < vEff_med_nd(bb)));



    % --- with drizzle ---
    re_med_d(bb)     = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)    = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)     = median(vertSeg_d{bb,3}, 'omitnan');
    vEff_med_d(bb)   = median(vertSeg_d{bb,4}, 'omitnan');

    re = vertSeg_d{bb,1};
    re_avgDevAbove_d(bb)    = mean( re(re > re_med_d(bb)) - re_med_d(bb));
    re_avgDevBelow_d(bb)    = mean( re_med_d(bb) - re(re < re_med_d(bb)));

    lwc = vertSeg_d{bb,2};
    lwc_avgDevAbove_d(bb)    = mean( lwc(lwc > lwc_med_d(bb)) - lwc_med_d(bb));
    lwc_avgDevBelow_d(bb)    = mean( lwc_med_d(bb) - lwc(lwc < lwc_med_d(bb)));

    Nc = vertSeg_d{bb,3};
    Nc_avgDevAbove_d(bb)    = mean( Nc(Nc > Nc_med_d(bb)) - Nc_med_d(bb));
    Nc_avgDevBelow_d(bb)    = mean( Nc_med_d(bb) - Nc(Nc < Nc_med_d(bb)));

    vEff = vertSeg_d{bb,4};
    vEff_avgDevAbove_d(bb)    = mean( vEff(vEff > vEff_med_d(bb)) - vEff_med_d(bb));
    vEff_avgDevBelow_d(bb)    = mean( vEff_med_d(bb) - vEff(vEff < vEff_med_d(bb)));

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

subplot(2,4,1)
hold on
x = [re_med_nd - re_avgDevBelow_nd; flipud(re_med_nd + re_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_avgDevBelow_d; flipud(re_med_d + re_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 17])

% Create textbox
annotation(figure1,'textbox',[0.130138 0.836 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');






subplot(2,4,2)
hold on
x = [vEff_med_nd - vEff_avgDevBelow_nd; flipud(vEff_med_nd + vEff_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [vEff_med_d - vEff_avgDevBelow_d; flipud(vEff_med_d + vEff_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.01, 0.175])

% Create textbox
annotation(figure1,'textbox',[0.451315023418542 0.836 0.0373 0.0666],'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');








subplot(2,4,3)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_avgDevBelow_nd; flipud(lwc_med_nd + lwc_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_avgDevBelow_d; flipud(lwc_med_d + lwc_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
% Create textbox
annotation(figure1,'textbox',[0.540879493143727 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 0.675])






subplot(2,4,4)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_avgDevBelow_nd; flipud(Nc_med_nd + Nc_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_avgDevBelow_d; flipud(Nc_med_d + Nc_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% Create textbox
annotation(figure1,'textbox',[0.746454846227638 0.836 0.0373 0.0666],'String','(d)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([0, 270])

legend({'','Non-drizzling',...
    '','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])
% legend({'',['Non-drizzling - N = ', num2str(length(idx_no_drizzle))],...
%     '',['Drizzling - N = ', num2str(length(idx_drizzle))]}, 'Location', 'best',...
%     'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
%     'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])




% Set the height and width as multiples of the IEEE paper requirements
w = 7.16; % inches
h = 3;    % inches

figure1.Units = 'inches';

figure1.Position = [1, 1, 2.75*w, 2.75*h]; % [left, bottom, width, height]










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
vertSeg_nd = cell(n_bins, 4);   % {re, lwc, Nc, var_eff}
vertSeg_d  = cell(n_bins, 4);

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


    % ------------------------------------------------------------------
    % Fit the droplet size distribution to get the gamma shape parameter
    % This is the ORACLES equivalent of the pre-stored gammaFit.alpha from
    % VOCALS-REx profiles.  We compute it here because find_verticalProfiles_ORACLES
    % does not store the fit result on the profile struct.

    significance_lvl = 0.1;

    [~, ~, gammaFit] = find_bestFitDist_dropDist( ...
        ensemble_profiles{nn}.Nd, ...
        ensemble_profiles{nn}.drop_radius_bin_edges, ...
        ensemble_profiles{nn}.drop_radius_bin_center, ...
        significance_lvl);

    v_eff_bin = 1 ./ (3 + gammaFit.alpha);            % per-altitude effective variance
    % -------------------------------------------------------------------------


    % Orient profile so index 1 is cloud top (tau = 0)
    % dz_dt > 0 → ascending (starts at base, tau increases going up)
    dz_dt = mean(diff(ensemble_profiles{nn}.altitude) ./ ...
        diff(ensemble_profiles{nn}.time));
    if dz_dt > 0
        % Ascending: flip so cloud top is first
        re_bin  = fliplr(ensemble_profiles{nn}.re);
        lwc_bin = fliplr(ensemble_profiles{nn}.lwc);
        Nc_bin  = fliplr(ensemble_profiles{nn}.total_Nc);
        v_eff_bin = flipud(v_eff_bin);
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
            vertSeg_d{bb,4} = [vertSeg_d{bb,4},  v_eff_bin(idx_seg)'];
        else
            vertSeg_nd{bb,1} = [vertSeg_nd{bb,1}, re_bin(idx_seg)];
            vertSeg_nd{bb,2} = [vertSeg_nd{bb,2}, lwc_bin(idx_seg)];
            vertSeg_nd{bb,3} = [vertSeg_nd{bb,3}, Nc_bin(idx_seg)];
            vertSeg_nd{bb,4} = [vertSeg_nd{bb,4}, v_eff_bin(idx_seg)'];
        end

    end

end




% Compute median and asymmetric average deviation (above/below median) for each bin
re_med_nd    = zeros(n_bins,1);  re_avgDevAbove_nd    = zeros(n_bins,1); re_avgDevBelow_nd    = zeros(n_bins,1);
lwc_med_nd   = zeros(n_bins,1);  lwc_avgDevAbove_nd   = zeros(n_bins,1); lwc_avgDevBelow_nd   = zeros(n_bins,1);
Nc_med_nd    = zeros(n_bins,1);  Nc_avgDevAbove_nd    = zeros(n_bins,1); Nc_avgDevBelow_nd    = zeros(n_bins,1);
vEff_med_nd  = zeros(n_bins,1);  vEff_avgDevAbove_nd  = zeros(n_bins,1); vEff_avgDevBelow_nd  = zeros(n_bins,1);

re_med_d     = zeros(n_bins,1);  re_avgDevAbove_d    = zeros(n_bins,1); re_avgDevBelow_d    = zeros(n_bins,1);
lwc_med_d    = zeros(n_bins,1);  lwc_avgDevAbove_d   = zeros(n_bins,1); lwc_avgDevBelow_d   = zeros(n_bins,1);
Nc_med_d     = zeros(n_bins,1);  Nc_avgDevAbove_d    = zeros(n_bins,1); Nc_avgDevBelow_d    = zeros(n_bins,1);
vEff_med_d   = zeros(n_bins,1);  vEff_avgDevAbove_d  = zeros(n_bins,1); vEff_avgDevBelow_d  = zeros(n_bins,1);

for bb = 1:n_bins
    % --- without drizzle ---
    re_med_nd(bb)    = median(vertSeg_nd{bb,1}, 'omitnan');
    lwc_med_nd(bb)   = median(vertSeg_nd{bb,2}, 'omitnan');
    Nc_med_nd(bb)    = median(vertSeg_nd{bb,3}, 'omitnan');
    vEff_med_nd(bb)  = median(vertSeg_nd{bb,4}, 'omitnan');

    re = vertSeg_nd{bb,1};
    re_avgDevAbove_nd(bb)    = mean( re(re > re_med_nd(bb)) - re_med_nd(bb));
    re_avgDevBelow_nd(bb)    = mean( re_med_nd(bb) - re(re < re_med_nd(bb)));

    lwc = vertSeg_nd{bb,2};
    lwc_avgDevAbove_nd(bb)    = mean( lwc(lwc > lwc_med_nd(bb)) - lwc_med_nd(bb));
    lwc_avgDevBelow_nd(bb)    = mean( lwc_med_nd(bb) - lwc(lwc < lwc_med_nd(bb)));

    Nc = vertSeg_nd{bb,3};
    Nc_avgDevAbove_nd(bb)    = mean( Nc(Nc > Nc_med_nd(bb)) - Nc_med_nd(bb));
    Nc_avgDevBelow_nd(bb)    = mean( Nc_med_nd(bb) - Nc(Nc < Nc_med_nd(bb)));

    vEff = vertSeg_nd{bb,4};
    vEff_avgDevAbove_nd(bb)    = mean( vEff(vEff > vEff_med_nd(bb)) - vEff_med_nd(bb));
    vEff_avgDevBelow_nd(bb)    = mean( vEff_med_nd(bb) - vEff(vEff < vEff_med_nd(bb)));


    % --- with drizzle ---
    re_med_d(bb)     = median(vertSeg_d{bb,1}, 'omitnan');
    lwc_med_d(bb)    = median(vertSeg_d{bb,2}, 'omitnan');
    Nc_med_d(bb)     = median(vertSeg_d{bb,3}, 'omitnan');
    vEff_med_d(bb)   = median(vertSeg_d{bb,4}, 'omitnan');

    re = vertSeg_d{bb,1};
    re_avgDevAbove_d(bb)    = mean( re(re > re_med_d(bb)) - re_med_d(bb));
    re_avgDevBelow_d(bb)    = mean( re_med_d(bb) - re(re < re_med_d(bb)));

    lwc = vertSeg_d{bb,2};
    lwc_avgDevAbove_d(bb)    = mean( lwc(lwc > lwc_med_d(bb)) - lwc_med_d(bb));
    lwc_avgDevBelow_d(bb)    = mean( lwc_med_d(bb) - lwc(lwc < lwc_med_d(bb)));

    Nc = vertSeg_d{bb,3};
    Nc_avgDevAbove_d(bb)    = mean( Nc(Nc > Nc_med_d(bb)) - Nc_med_d(bb));
    Nc_avgDevBelow_d(bb)    = mean( Nc_med_d(bb) - Nc(Nc < Nc_med_d(bb)));

    vEff = vertSeg_d{bb,4};
    vEff_avgDevAbove_d(bb)    = mean( vEff(vEff > vEff_med_d(bb)) - vEff_med_d(bb));
    vEff_avgDevBelow_d(bb)    = mean( vEff_med_d(bb) - vEff(vEff < vEff_med_d(bb)));

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

subplot(2,4,5)
hold on
x = [re_med_nd - re_avgDevBelow_nd; flipud(re_med_nd + re_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [re_med_d - re_avgDevBelow_d; flipud(re_med_d + re_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(re_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 17])

% Create textbox
annotation(figure1,'textbox',[0.132465746741153 0.367111111111112 0.0373 0.0666],'String',{'(e)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');





subplot(2,4,6)
hold on
x = [vEff_med_nd - vEff_avgDevBelow_nd; flipud(vEff_med_nd + vEff_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [vEff_med_d - vEff_avgDevBelow_d; flipud(vEff_med_d + vEff_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(vEff_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
xlabel('$\langle \nu_{e}(\tau) \rangle$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.01, 0.175])

% Create textbox
annotation(figure1,'textbox',[0.343867473506009 0.367111111111112 0.0373 0.0666],'String',{'(f)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');









subplot(2,4,7)

% axes2 = axes('Parent',figure1,'Position',[0.465797101449277 0.11 0.21340579710145 0.815]);
% hold(axes2,'on');

hold on
x = [lwc_med_nd - lwc_avgDevBelow_nd; flipud(lwc_med_nd + lwc_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [lwc_med_d - lwc_avgDevBelow_d; flipud(lwc_med_d + lwc_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(lwc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.675])

% Create textbox
annotation(figure1,'textbox',[0.542078635404322 0.361555555555556 0.0372999999999999 0.0666],'String','(g)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');







subplot(2,4,8)

% Create axes
% axes3 = axes('Parent',figure1,'Position',[0.746594202898554 0.11 0.213405797101449 0.815]);
% hold(axes3,'on');

hold on
x = [Nc_med_nd - Nc_avgDevBelow_nd; flipud(Nc_med_nd + Nc_avgDevAbove_nd)];
y = [bin_center; flipud(bin_center)];
fill(x, y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_nd, bin_center, '-', 'Color', plt_clr_1, 'LineWidth', 1.5)
x = [Nc_med_d - Nc_avgDevBelow_d; flipud(Nc_med_d + Nc_avgDevAbove_d)];
fill(x, y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', 0.2)
plot(Nc_med_d, bin_center, '-', 'Color', plt_clr_2, 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
% ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 270])

% Create textbox
annotation(figure1,'textbox',[0.745932819705434 0.367111111111112 0.0373 0.0666],'String','(h)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% legend({'',['Non-drizzling - N = ', num2str(length(idx_no_drizzle))],...
%     '',['Drizzling - N = ', num2str(length(idx_drizzle))]}, 'Location', 'best',...
%     'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
%     'Position',[0.871976270152657 0.367111111111112 0.122791322836865 0.111660079616977])





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
    '- 4 panels with logNorm uncert.fig']);


% save .png with 500 DPI resolution
% remove title
title('');
exportgraphics(figure1,[folderpath_figs,'VOCALS-REx and ORACLES in-situ median profiles separated by drizzle',...
    'and non drizzle - 4 panels with logNorm uncert.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------








%% Combine VOCALS-REx and ORACLES - percentile bands (10-90 outer, 25-75 inner)
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
%   outer (lighter) band : 10th to 90th percentile
%   inner (darker)  band : 25th to 75th percentile
%   solid line           : median (50th percentile)
%
% Because each edge of each band is the actual percentile of the sampled
% data at that vertical level, the left/right widths relative to the
% median reflect the true asymmetry in the distribution.


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

    tau_prof = ensemble_profiles{nn}.tau;
    if max(tau_prof) < tauC_limit
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

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
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
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
pct = [10 25 50 75 90];
re_p_nd   = zeros(n_bins, 5);  re_p_d   = zeros(n_bins, 5);
lwc_p_nd  = zeros(n_bins, 5);  lwc_p_d  = zeros(n_bins, 5);
Nc_p_nd   = zeros(n_bins, 5);  Nc_p_d   = zeros(n_bins, 5);
vEff_p_nd = zeros(n_bins, 5);  vEff_p_d = zeros(n_bins, 5);

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

plt_clr_1   = mySavedColors(64,'fixed');
plt_clr_2   = mySavedColors(62,'fixed');

fnt_sz      = 18;
ttl_fnt     = 20;
alpha_outer = 0.17;
alpha_inner = 0.30;
ln_w        = 1.5;


y = [bin_center; flipud(bin_center)];

% --- (a) r_e ---
subplot(2,4,1)
hold on
fill([re_p_nd(:,1); flipud(re_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([re_p_nd(:,2); flipud(re_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(re_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([re_p_d(:,1);  flipud(re_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([re_p_d(:,2);  flipud(re_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(re_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
grid on; grid minor
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 18])

annotation(figure1,'textbox',[0.130138 0.836 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (b) v_eff ---
subplot(2,4,2)
hold on
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([vEff_p_nd(:,2); flipud(vEff_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(vEff_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([vEff_p_d(:,2);  flipud(vEff_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(vEff_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
xlim([0.005, 0.175])

annotation(figure1,'textbox',[0.451315023418542 0.836 0.0373 0.0666],'String',{'(b)'},...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (c) LWC ---
subplot(2,4,3)
hold on
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([lwc_p_nd(:,2); flipud(lwc_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(lwc_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([lwc_p_d(:,2);  flipud(lwc_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(lwc_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
grid on; grid minor
title(['VOCALS-REx median profiles - N = ', num2str(num_kept),' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
annotation(figure1,'textbox',[0.540879493143727 0.836 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');
xlim([0, 0.7])

% --- (d) N_c ---
subplot(2,4,4)
hold on
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([Nc_p_nd(:,2); flipud(Nc_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
h_nd = plot(Nc_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w);
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([Nc_p_d(:,2);  flipud(Nc_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
h_d  = plot(Nc_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w);
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlim([0, 300])
annotation(figure1,'textbox',[0.746454846227638 0.836 0.0373 0.0666],'String','(d)',...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

legend([h_nd h_d], {'Non-drizzling','Drizzling'}, ...
    'Location','best','Interpreter','latex', ...
    'FontSize', fnt_sz-2, 'Color','white','TextColor','black', ...
    'Position',[0.871976270152657 0.82592592592596 0.122791322836865 0.111660079616977])


% -------------------------------------------------------------------------
% ------ ORACLES: load, separate drizzle, compute percentiles -------------
% -------------------------------------------------------------------------
clearvars -except figure1 bin_center n_bins bin_edges fnt_sz ttl_fnt ...
    plt_clr_1 plt_clr_2 alpha_outer alpha_inner ln_w

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

    tau_prof = ensemble_profiles{nn}.tau;
    if max(tau_prof) < tauC_limit
        num_removed = num_removed + 1;
        continue
    else
        num_kept = num_kept + 1;
    end
    tau_norm = tau_prof ./ max(tau_prof);

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
            idx_seg = tau_norm >= bin_edges(bb) & tau_norm <= bin_edges(bb+1);
        else
            idx_seg = tau_norm >  bin_edges(bb) & tau_norm <= bin_edges(bb+1);
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

pct = [10 25 50 75 90];
re_p_nd   = zeros(n_bins, 5);  re_p_d   = zeros(n_bins, 5);
lwc_p_nd  = zeros(n_bins, 5);  lwc_p_d  = zeros(n_bins, 5);
Nc_p_nd   = zeros(n_bins, 5);  Nc_p_d   = zeros(n_bins, 5);
vEff_p_nd = zeros(n_bins, 5);  vEff_p_d = zeros(n_bins, 5);

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
fill([re_p_nd(:,1); flipud(re_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([re_p_nd(:,2); flipud(re_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(re_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([re_p_d(:,1);  flipud(re_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([re_p_d(:,2);  flipud(re_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(re_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle r_e(\tau) \rangle \; (\mu m)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
ylabel('Normalized Optical Depth', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([3, 18])
annotation(figure1,'textbox',[0.132465746741153 0.367111111111112 0.0373 0.0666],'String',{'(e)'},...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (f) v_eff ---
subplot(2,4,6)
hold on
fill([vEff_p_nd(:,1); flipud(vEff_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([vEff_p_nd(:,2); flipud(vEff_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(vEff_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([vEff_p_d(:,1);  flipud(vEff_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([vEff_p_d(:,2);  flipud(vEff_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(vEff_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
set(gca, 'XScale', 'log')
grid on; grid minor
xlabel('$\langle \nu_{e}(\tau) \rangle$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0.005, 0.175])
annotation(figure1,'textbox',[0.343867473506009 0.367111111111112 0.0373 0.0666],'String',{'(f)'},...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (g) LWC ---
subplot(2,4,7)
hold on
fill([lwc_p_nd(:,1); flipud(lwc_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([lwc_p_nd(:,2); flipud(lwc_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(lwc_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([lwc_p_d(:,1);  flipud(lwc_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([lwc_p_d(:,2);  flipud(lwc_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(lwc_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle LWC(\tau) \rangle \; (g/m^3)$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
title(['ORACLES median profiles - N = ', num2str(num_kept), ' - $\tau_{c} \geq$ ', num2str(tauC_limit)],...
    'Interpreter', 'latex','FontSize', ttl_fnt)
xlim([0, 0.7])
annotation(figure1,'textbox',[0.542078635404322 0.361555555555556 0.0372999999999999 0.0666],'String','(g)',...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% --- (h) N_c ---
subplot(2,4,8)
hold on
fill([Nc_p_nd(:,1); flipud(Nc_p_nd(:,5))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([Nc_p_nd(:,2); flipud(Nc_p_nd(:,4))], y, plt_clr_1, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(Nc_p_nd(:,3), bin_center, '-', 'Color', plt_clr_1, 'LineWidth', ln_w)
fill([Nc_p_d(:,1);  flipud(Nc_p_d(:,5))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_outer)
fill([Nc_p_d(:,2);  flipud(Nc_p_d(:,4))],  y, plt_clr_2, 'EdgeAlpha', 0, 'FaceAlpha', alpha_inner)
plot(Nc_p_d(:,3),  bin_center, '-', 'Color', plt_clr_2, 'LineWidth', ln_w)
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$\langle N_c(\tau) \rangle \; (cm^{-3})$', 'Interpreter', 'latex', 'FontSize', fnt_sz)
xlim([0, 300])
annotation(figure1,'textbox',[0.745932819705434 0.367111111111112 0.0373 0.0666],'String','(h)',...
    'Interpreter','latex','FontSize',25,'FontName','Helvetica Neue',...
    'FitBoxToText','off','EdgeColor','none');

% IEEE-column sizing
w = 7.16; h = 3;
figure1.Units    = 'inches';
figure1.Position = [1, 1, 2.75*w, 2.75*h];




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











%% Create imagesc of the droplet distribution from VOCALS-REx and ORACLES with overlaid median profile

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

