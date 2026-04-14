


%% Create median profiles of re, lwc and tau_c from ORACLES data

clear variables

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

% Search for ensemble profiles mat file (use the most recently created one)
folder_contents = dir(foldername_data);
mat_candidates = {};

for nn = 1:length(folder_contents)
    fname = folder_contents(nn).name;
    if length(fname) > 17 && strcmp(fname(1:17), 'ensemble_profiles') && ...
            strcmp(fname(end-3:end), '.mat')
        mat_candidates{end+1} = fname; %#ok<SAGROW>
    end
end

if isempty(mat_candidates)
    error('No ensemble_profiles mat file found in:\n  %s', foldername_data)
end

% Use the last file listed (typically most recent by name/date)
mat_filename = mat_candidates{end};
disp(['Loading: ', mat_filename])
load([foldername_data, mat_filename], 'ensemble_profiles', 'inputs')

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
annotation(figure1,'textbox',[0.191 0.806 0.0373 0.0666],'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.00447158452278369 0.806666666666669 0.14136174881055 0.111660079616977])







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
title('ORACLES median profiles', 'Interpreter', 'latex',...
    'FontSize', ttl_fnt)
% Create textbox
annotation(figure1,'textbox',[0.469 0.806 0.0373 0.0666],'String','(b)',...
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
annotation(figure1,'textbox',[0.745 0.806 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

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
% saveas(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(figure1,[folderpath_figs,'ORACLES in-situ profiles separated by drizzle and non drizzle',...
%     '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------







%% %% Create median profiles of re, lwc and tau_c from VOCALS-REx data

clear variables

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
vertSeg_nd = cell(n_bins, 3);   % {re, lwc, Nc}
vertSeg_d  = cell(n_bins, 3);

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

legend({'','Non-drizzling','','Drizzling'}, 'Location', 'best',...
    'Interpreter', 'latex', 'FontSize', fnt_sz-2, 'Color', 'white', 'TextColor', 'black',...
    'Position',[0.00447158452278369 0.806666666666669 0.14136174881055 0.111660079616977])

% Create textbox
annotation(figure1,'textbox',[0.191 0.806 0.0373 0.0666],'String',{'(a)'},...
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
annotation(figure1,'textbox',[0.469 0.806 0.0373 0.0666],'String','(b)',...
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
annotation(figure1,'textbox',[0.745 0.806 0.0373 0.0666],'String','(c)',...
    'Interpreter','latex',...
    'FontSize',25,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

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

