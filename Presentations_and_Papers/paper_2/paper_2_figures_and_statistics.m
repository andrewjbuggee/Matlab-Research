%% Figures for Paper 2

% By Andrew John Buggee




%% Figure: Impact of IWV assumption on droplet size retrieval bias
%
% For each of ~840 simulated cloud retrievals, compute the signed retrieval
% bias introduced by each of four assumed total IWV values (10, 15, 20, 25 mm)
% relative to the true (simulated) cloud state:
%
%   delta_r_top(i,j) = r_top_retrieved(i,j) - r_top_true(i)
%   delta_r_bot(i,j) = r_bot_retrieved(i,j) - r_bot_true(i)
%
% Panel (a): Scatter of delta_r_top vs. true IWV_ac (above-cloud PW)
%            with binned medians +/- IQR. One color per IWV assumption.
% Panel (b): Same for delta_r_bot vs. true IWV_ac.
% Panel (c): Box plots of the sensitivity d(delta_r)/d(IWV_assumed)
%            (slope of bias vs. assumed IWV fitted across 4 assumptions),
%            separately for r_top and r_bot, binned by cloud optical depth.

clear variables

% -------------------------------------------------------------------------
% Setup: folder paths for each assumed IWV retrieval
% -------------------------------------------------------------------------
% which_computer = whatComputer();

base_path = ['/Users/', whatComputer, '/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    'paper2_variableSweep/dropProf_retrieval_logTransform_newCov_noACPW_retrieval_assumed_10_15_20_25'];


folder_names = { ...
    'newRetrieval_logTransform_newCov_noACPW_10', ...
    'newRetrieval_logTransform_newCov_noACPW_15', ...
    'newRetrieval_logTransform_newCov_noACPW_20_rev2', ...
    'newRetrieval_logTransform_newCov_noACPW_25_rev2' };

iwv_totalColumn_assumed = [10, 15, 20, 25];   % assumed total column water vapor (mm)
iwv_aboveCloud_assumed = [5.7864, 8.6797, 11.5729, 14.4661];
n_iwv       = numel(iwv_aboveCloud_assumed);

% -------------------------------------------------------------------------
% Load GN retrieval data from all four folders
% Only _1.mat files (first noise realization) with a GN_outputs struct are
% used; duplicate trial runs and incomplete files are silently skipped.
% -------------------------------------------------------------------------
fprintf('Loading retrieval data...\n')

rtop_true_all  = cell(1, n_iwv);
rbot_true_all  = cell(1, n_iwv);
tauc_true_all  = cell(1, n_iwv);
actpw_true_all = cell(1, n_iwv);
rtop_ret_all   = cell(1, n_iwv);
rbot_ret_all   = cell(1, n_iwv);

for jj = 1:n_iwv

    folder = fullfile(base_path, folder_names{jj});
    files  = dir(fullfile(folder, 'dropletRetrieval_noACPW_*_1.mat'));

    % Pre-allocate with maximum possible size
    n         = 0;
    N_max     = numel(files);
    rtop_t    = nan(N_max, 1);
    rbot_t    = nan(N_max, 1);
    tauc_t    = nan(N_max, 1);
    actpw_t   = nan(N_max, 1);
    rtop_r    = nan(N_max, 1);
    rbot_r    = nan(N_max, 1);

    for nn = 1:N_max
        fpath = fullfile(files(nn).folder, files(nn).name);
        ds    = load(fpath);

        if ~isfield(ds, 'GN_outputs'), continue; end

        n = n + 1;
        rtop_t(n)  = ds.GN_inputs.measurement.r_top;
        rbot_t(n)  = ds.GN_inputs.measurement.r_bot;
        tauc_t(n)  = ds.GN_inputs.measurement.tau_c;
        actpw_t(n) = ds.GN_inputs.measurement.actpw;
        rtop_r(n)  = ds.GN_outputs.retrieval(1, end);  % row 1 = r_top
        rbot_r(n)  = ds.GN_outputs.retrieval(2, end);  % row 2 = r_bot
    end

    rtop_true_all{jj}  = rtop_t(1:n);
    rbot_true_all{jj}  = rbot_t(1:n);
    tauc_true_all{jj}  = tauc_t(1:n);
    actpw_true_all{jj} = actpw_t(1:n);
    rtop_ret_all{jj}   = rtop_r(1:n);
    rbot_ret_all{jj}   = rbot_r(1:n);

    fprintf('  IWV = %f mm : %d retrievals loaded\n', iwv_aboveCloud_assumed(jj), n)

end

% Signed retrieval biases
delta_rtop_all = cellfun(@(r,t) r - t, rtop_ret_all, rtop_true_all, 'UniformOutput', false);
delta_rbot_all = cellfun(@(r,t) r - t, rbot_ret_all, rbot_true_all, 'UniformOutput', false);

% -------------------------------------------------------------------------
% Figure layout and style
% -------------------------------------------------------------------------
C            = mySavedColors([66:69], 'fixed');   % one color per IWV assumption
fnt_sz       = 18;
dot_sz       = 18;
alpha_sc     = 0.20;    % scatter point transparency
lw           = 2.5;     % line width for median curves
mk_sz        = 30;

% x-axis bins for actpw (above-cloud precipitable water)
% actpw_edges = 3 : 1 : 12;          % 1 mm bins over [3, 12] mm range
% actpw_edges = [3, 4.25, 5.5, 6.75, 8, 9.25, 11];          % custom bins to fit the true value modled in this data set
% actpw_ctrs  = actpw_edges(1:end-1) + 0.5;
actpw_ctrs = unique(actpw_true_all{1})';                     % The true centers of the 6 different ACPWs modeled in this data set
actpw_edges = [actpw_ctrs - 0.25, actpw_ctrs(end) + 0.25];          % custom bins to fit the true value modled in this data set
n_bins_pw   = numel(actpw_ctrs);

fig1 = figure('Position', [50 50 1600 520]);

% -------------------------------------------------------------------------
% Panels (a) and (b): scatter + binned medians with IQR shading
% -------------------------------------------------------------------------
panel_titles    = {'(a) Cloud-top radius bias', '(b) Cloud-base radius bias'};
delta_list      = {delta_rtop_all, delta_rbot_all};
ylabels_list    = {'$\Delta r_{top}$ ($\mu$m)', '$\Delta r_{bot}$ ($\mu$m)'};

for pp = 1:2

    subplot(1, 2, pp);  hold on;  box on;  grid on;  grid minor

    delta_all = delta_list{pp};

    % --- raw scatter (context) ---
    % for jj = 1:n_iwv
    %     scatter(actpw_true_all{jj}, delta_all{jj}, dot_sz, ...
    %         'MarkerFaceColor', C(jj,:), 'MarkerEdgeColor', 'none', ...
    %         'MarkerFaceAlpha', alpha_sc)
    % end

    % --- binned medians + IQR shading ---
    h_lines = gobjects(n_iwv, 1);
    for jj = 1:n_iwv
        med_bias = nan(1, n_bins_pw);
        q25      = nan(1, n_bins_pw);
        q75      = nan(1, n_bins_pw);

        for bb = 1:n_bins_pw
            in_bin = actpw_true_all{jj} >= actpw_edges(bb) & ...
                     actpw_true_all{jj} <  actpw_edges(bb+1);
            if sum(in_bin) > 0
                med_bias(bb) = median(delta_all{jj}(in_bin));
                q25(bb)      = prctile(delta_all{jj}(in_bin), 25);
                q75(bb)      = prctile(delta_all{jj}(in_bin), 75);
            end
        end

        valid = ~isnan(med_bias);

        fill([actpw_ctrs(valid), fliplr(actpw_ctrs(valid))], ...
             [q25(valid),        fliplr(q75(valid))], ...
             C(jj,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off')

        h_lines(jj) = plot(actpw_ctrs(valid), med_bias(valid), '.-', ...
            'Color', C(jj,:), 'LineWidth', lw, 'MarkerSize', mk_sz,...
            'DisplayName', sprintf('Assumed IWV = %2.1f mm', iwv_aboveCloud_assumed(jj)));

        % h_lines(jj) = errorbar(actpw_ctrs(valid), med_bias(valid), q25(valid), q75(valid),  '.-', ...
        %     'Color', C(jj,:), 'LineWidth', lw, 'MarkerSize', mk_sz,...
        %     'DisplayName', sprintf('Assumed IWV = %2.1f mm', iwv_aboveCloud_assumed(jj)));
    end

    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')

    xlabel('True $IWV_{ac}$ (mm)', 'Interpreter', 'latex', 'FontSize', fnt_sz)
    ylabel(ylabels_list{pp},       'Interpreter', 'latex', 'FontSize', fnt_sz)
    title(panel_titles{pp},        'Interpreter', 'latex', 'FontSize', fnt_sz)

    if pp == 1
        legend(h_lines, 'Interpreter', 'latex', 'FontSize', 19, 'Location', 'best')
    end

    set(gca, 'FontSize', fnt_sz)

end

% -------------------------------------------------------------------------
% Overall title and final sizing
% -------------------------------------------------------------------------
sgtitle('Impact of IWV assumption on droplet size retrieval', ...
    'Interpreter', 'latex', 'FontSize', 16)

set(gcf, 'Color', 'w')

% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig1,[folderpath_figs,'Retrieval bias due to assumed IWV content above cloud.fig']);


% save .png with 500 DPI resolution
% remove title
subplot(1,2,1); title(''); subplot(1,2,2); title(''); sgtitle('');
exportgraphics(fig1,[folderpath_figs,'Retrieval bias due to assumed IWV content above cloud.jpg'],...
    'Resolution', 500);
% -------------------------------------
% -------------------------------------





% -------------------------------------------------------------------------
% Panel (c): sensitivity d(delta_r)/d(IWV_assumed), binned by tau_c
%
% For each cloud present in all 4 assumption folders, fit a line through the
% 4 (IWV_assumed, delta_r) points and record the slope as the sensitivity.
% Use IWV=20 folder (index 3, 840 complete retrievals) as the reference set.
% -------------------------------------------------------------------------
ref = 3;        % index of the reference folder (IWV=20)
N_ref = numel(actpw_true_all{ref});
tol   = 5e-3;   % numerical tolerance for matching truth values

sens_rtop = nan(N_ref, 1);
sens_rbot = nan(N_ref, 1);
tauc_sens = nan(N_ref, 1);

for ii = 1:N_ref

    rt = rtop_true_all{ref}(ii);
    rb = rbot_true_all{ref}(ii);
    tc = tauc_true_all{ref}(ii);
    ap = actpw_true_all{ref}(ii);

    d_top  = nan(1, n_iwv);
    d_bot  = nan(1, n_iwv);

    for jj = 1:n_iwv
        idx = find( abs(rtop_true_all{jj} - rt) < tol & ...
                    abs(rbot_true_all{jj} - rb) < tol & ...
                    abs(tauc_true_all{jj} - tc) < tol & ...
                    abs(actpw_true_all{jj} - ap) < tol, 1);
        if ~isempty(idx)
            d_top(jj) = delta_rtop_all{jj}(idx);
            d_bot(jj) = delta_rbot_all{jj}(idx);
        end
    end

    % Require at least 3 of 4 matched assumptions to fit the slope
    valid_mask = ~isnan(d_top);
    if sum(valid_mask) >= 3
        p_top = polyfit(double(iwv_aboveCloud_assumed(valid_mask)), d_top(valid_mask), 1);
        p_bot = polyfit(double(iwv_aboveCloud_assumed(valid_mask)), d_bot(valid_mask), 1);
        sens_rtop(ii) = p_top(1);   % units: um / mm
        sens_rbot(ii) = p_bot(1);
        tauc_sens(ii) = tc;
    end

end

fprintf('Sensitivity computed for %d matched clouds\n', sum(~isnan(sens_rtop)))

% tau_c bins
tauc_edges       = [6, 12, 18, 25];
n_tauc_bins      = numel(tauc_edges) - 1;
tauc_bin_labels  = {'6 \leq \tau_c < 12', '12 \leq \tau_c < 18', '18 \leq \tau_c \leq 24'};

% Colors for r_top and r_bot boxes
c_top = [0.22, 0.45, 0.77];   % blue
c_bot = [0.80, 0.25, 0.18];   % red

% Build a padded matrix (NaN-filled) for boxplot:
%   odd columns  (1, 3, 5) = r_top sensitivity in each tau bin
%   even columns (2, 4, 6) = r_bot sensitivity in each tau bin
box_data_top = cell(1, n_tauc_bins);
box_data_bot = cell(1, n_tauc_bins);

for bb = 1:n_tauc_bins
    in_bin = tauc_sens >= tauc_edges(bb) & tauc_sens < tauc_edges(bb+1);
    box_data_top{bb} = sens_rtop(in_bin & ~isnan(sens_rtop));
    box_data_bot{bb} = sens_rbot(in_bin & ~isnan(sens_rbot));
end

max_n = max(cellfun(@numel, [box_data_top, box_data_bot]));
bmat  = nan(max_n, 2 * n_tauc_bins);
for bb = 1:n_tauc_bins
    bmat(1:numel(box_data_top{bb}), 2*bb-1) = box_data_top{bb};
    bmat(1:numel(box_data_bot{bb}), 2*bb  ) = box_data_bot{bb};
end

% Box positions: pair of boxes per bin, bins spread 3 units apart
bin_center = (1:n_tauc_bins) * 3;
pos = zeros(1, 2 * n_tauc_bins);
for bb = 1:n_tauc_bins
    pos(2*bb-1) = bin_center(bb) - 0.5;   % r_top
    pos(2*bb  ) = bin_center(bb) + 0.5;   % r_bot
end

subplot(1, 3, 3);  hold on;  box on;  grid on

boxplot(bmat, 'Positions', pos, 'Widths', 0.75, 'Symbol', 'k+', ...
    'Colors', repmat([0 0 0], 2*n_tauc_bins, 1))

% Recolor box fill: findobj returns handles in reverse draw order
% (last drawn = first element), so h(1) = rightmost box = r_bot_bin3
h_boxes = findobj(gca, 'Tag', 'Box');
n_boxes = numel(h_boxes);
for kk = 1:n_boxes
    group_idx = n_boxes - kk + 1;  % map from reverse order to forward (1..6)
    if mod(group_idx, 2) == 1      % odd columns = r_top
        fc = c_top;
    else                           % even columns = r_bot
        fc = c_bot;
    end
    patch(get(h_boxes(kk), 'XData'), get(h_boxes(kk), 'YData'), ...
        fc, 'FaceAlpha', 0.70, 'EdgeColor', 'k', 'LineWidth', 1.5)
end

yline(0, 'k--', 'LineWidth', 1.5)

% x-tick labels at bin centers
set(gca, 'XTick', bin_center, 'XTickLabel', tauc_bin_labels, ...
    'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz)

ylabel('$\partial(\Delta r) / \partial IWV_{\mathrm{assumed}}$ ($\mu$m mm$^{-1}$)', ...
    'Interpreter', 'latex', 'FontSize', fnt_sz)
title('(c) Retrieval sensitivity vs.\ $\tau_c$', 'Interpreter', 'latex', 'FontSize', fnt_sz)

% Legend patches
p_top_leg = patch(nan, nan, c_top, 'FaceAlpha', 0.70, 'EdgeColor', 'k', 'LineWidth', 1.5);
p_bot_leg = patch(nan, nan, c_bot, 'FaceAlpha', 0.70, 'EdgeColor', 'k', 'LineWidth', 1.5);
legend([p_top_leg, p_bot_leg], {'$r_{top}$', '$r_{bot}$'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')

set(gca, 'FontSize', fnt_sz)

% -------------------------------------------------------------------------
% Overall title and final sizing
% -------------------------------------------------------------------------
sgtitle('Impact of IWV assumption on droplet size retrieval', ...
    'Interpreter', 'latex', 'FontSize', 16)

set(gcf, 'Color', 'w')




%% Make paneled figure showing how the assumption of total water vapor column impacts the retrieval of a droplet profile

% *** 0.1% uncertainty ***

clear variables


% load 0.1% uncertainty data
filenames_noACPW_10 = 'newRetrieval_logTransform_newCov_noACPW_10';
filenames_noACPW_15 = 'newRetrieval_logTransform_newCov_noACPW_15';
filenames_noACPW_20 = 'newRetrieval_logTransform_newCov_noACPW_20';
filenames_noACPW_25 = 'newRetrieval_logTransform_newCov_noACPW_25';

folder_path = ['/Users/', whatComputer, '/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
    'HySICS/Droplet_profile_retrievals/vza_7_subset_pt1percent/'];



% what are the free parameters?
r_top = 10;
r_bot = 5;
% tau_c = [5,11,17,23];
tau_c = [5,11];
% tcpw = [8, 14, 20];
tcpw = [14];

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_tauC = length(tau_c);
num_tcpw = length(tcpw);




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');




for pw = 1:length(tcpw)

    f = figure;


    for tc = 1:length(tau_c)

        subplot(1,length(tau_c),tc)



        % Load all 5 retrievals
        fileNames = [dir([folder_path, filenames_noACPW_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);...
            dir([folder_path, filenames_full_startsWith, '_tauC_', num2str(tau_c(tc)),...
            '_tcwv_',  num2str(tcpw(pw)), '*.mat']);];

        lgnd_str = cell(1, length(fileNames)+2);

        for nn = 1:length(fileNames)



            % Load a data set
            ds = load([fileNames(nn).folder, '/', fileNames(nn).name]);


            if nn==1


                % first, plot the simulated profile values as two lines
                xline(ds.GN_inputs.measurement.r_bot, ':', ['Simulated $r_{bot}$'],...
                    'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255],...
                    'LabelHorizontalAlignment','left',...
                    'LabelVerticalAlignment','bottom');

                hold on

                yline(ds.GN_inputs.measurement.r_top, ':', ['Simulated $r_{top}$'],...
                    'Fontsize',20, 'Interpreter','latex','LineWidth',2,'Color', [147/255, 150/255, 151/255],...
                    'LabelHorizontalAlignment','right',...
                    'LabelVerticalAlignment','top');
                hold on




                % Skip the first two legend entries
                lgnd_str{1} = '';
                lgnd_str{2} = '';

                % write titles and subtitles
                if tc==1

                    % what was the assumed above cloud column water vapor path?
                    title(['Simulated measurements with $0.1\%$ uncertainty - $acpw$ = ',num2str(round(ds.GN_inputs.measurement.actpw, 2)), ' $mm$'],...
                        'Fontsize', 25, 'Interpreter', 'latex');

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.25 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==2


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.45 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==3


                    % Create textbox
                    annotation(f,'textbox',...
                        [0.66 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');

                elseif tc==4

                    % Create textbox
                    annotation(f,'textbox',...
                        [0.87 0.841666666666667 0.0369591836734694 0.0716666666666671],...
                        'String',{['$\tau_c = $ ', num2str(tau_c(tc))]},...
                        'Interpreter','latex',...
                        'FontSize',23,...
                        'FitBoxToText','off',...
                        'EdgeColor','none');



                end






            end



            % plot the retrieved droplet profile
            if nn<length(fileNames)

                hold on


                % Plot the retrieval uncertainty of the radius at cloud top and
                % bottom
                e1 = errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                    sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                    sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
                    'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                    'Color', C(nn+1,:));

                e1.Bar.LineStyle = 'dotted';

                hold on



            else

                % give a different marker type for the retrieval using 66 bands
                % that also retrieved above cloud column water vapor


                errorbar(ds.GN_outputs.retrieval(2,end), ds.GN_outputs.retrieval(1,end), sqrt(ds.GN_outputs.posterior_cov(1,1))/2,...
                    sqrt(ds.GN_outputs.posterior_cov(1,1))/2, sqrt(ds.GN_outputs.posterior_cov(2,2))/2,...
                    sqrt(ds.GN_outputs.posterior_cov(2,2))/2, 'MarkerFaceColor', C(nn+1,:),...
                    'MarkerEdgeColor', C(nn+1,:), 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 40,...
                    'Color', C(nn+1,:))

                hold on


            end



            % create the legend string
            if nn<length(fileNames)

                % what was the assumed above cloud column water vapor path?
                assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

                lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - assumed $acpw$ = ',...
                    num2str(round(assumed_CWV,2)), ' $mm$'];





            else

                % create the string for the retrieval using CWV

                % what was the retrieved above cloud column water vapor path above
                % cloud?
                retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

                lgnd_str{nn+2} = [num2str(numel(ds.GN_inputs.bands2run)), ' bands - retrieved $acpw$ = ',...
                    num2str(round(retrieved_CWV, 2)), ' $mm$'];

            end





        end






        % Create a Legend with only the two black curves
        legend(lgnd_str, 'Interpreter','latex', 'FontSize', 18,...
            'Position',[0.166173414792696 0.16735419938909 0.142261290258291 0.153166666030884])

        grid on; grid minor
        ylabel('$r_{top}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$r_{bot}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

        % set x and y lims
        % set figure limits
        dx = 2;
        dy = 1;

        xlim([ds.GN_inputs.measurement.r_bot - dx, ds.GN_inputs.measurement.r_bot + dx])
        ylim([ds.GN_inputs.measurement.r_top - dy, ds.GN_inputs.measurement.r_top + dy])



    end

    if strcmp(whatComputer, 'anbu8374')==true

        set(gcf,'Position',[0 0 2450 600])

    elseif strcmp(whatComputer, 'andrewbuggee')==true

        set(gcf,'Position',[0 0 1200 600])

    end

end









%% Plot the HySICS retrieval along with the in-situ measurement

% ** only considering re_profile uncertainty **

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];

    % define the folder where retrievals are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if strcmp(filenames_retrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];


% ------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_3 ---
% ------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
%plt_idx = 17;
% ------------------------------------------------------------


% -------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1 ---
% -------------------------------------------------------------------------------
% profile_indexes for paper = [3, 6, 7, 9, 18]
% plt_idx = 4;
% ------------------------------------------------------------


% ---------------------------------------------------------------------------------------------
% -- For full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3 ---
% ---------------------------------------------------------------------------------------------
% profile_indexes for paper = [2, 4, 10]
plt_idx = 4;
% ------------------------------------------------------------


fig1 = plot_retrieved_prof_with_inSitu_paper2(folder_paths.retrieval, filenames_retrieval(plt_idx).name);



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
% end
% saveas(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
%     num2str(plt_idx), '.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% title('');
% exportgraphics(fig1,[folderpath_figs,'HySICS retrieval with VR in-situ measurement - profile number ',...
%     num2str(plt_idx), '.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------



%% Compute HySICS retrieval statistics



clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    % % define the folder where retrievals are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_3/'];


    % define the folder where retrievals are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_and_cloudTop_uncert_3/'];

    addpath(folder_paths.retrieval);



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if strcmp(filenames_retrieval(nn).name(1), 'd')~=true

        idx_2delete = [idx_2delete, nn];

    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];



con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3



% store the LWP retrieval
% store the TBLUT LWP estimate
% store the TBLUT LWP estimate with Wood-Hartmann adjustement
% store the in-situ LWP measurement
lwp_retrieval = zeros(size(filenames_retrieval));
lwp_tblut = zeros(size(filenames_retrieval));
lwp_tblut_WH = zeros(size(filenames_retrieval));
lwp_inSitu = zeros(size(filenames_retrieval));

lwp_newCalc = zeros(size(filenames_retrieval));


% store the ACPW retrieval
% store the true ACPW used in the forward model - measured by....
acpw_retrieval = zeros(size(filenames_retrieval));
acpw_inSitu = zeros(size(filenames_retrieval));


% store the optical depth retrieval
% store the in-situ measured optical depth
tauC_retrieval = zeros(size(filenames_retrieval));
tauC_inSitu = zeros(size(filenames_retrieval));



for nn = 1:length(filenames_retrieval)

    clear ds

    ds = load(filenames_retrieval(nn).name);

    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the liquid water paths
    lwp_retrieval(nn) = ds.GN_outputs.LWP;    % g/m^2

    % compute the LWP estimate using the TBLUT retrieval
    lwp_tblut(nn) = (2 * rho_h2o * (ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

    % ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
    lwp_tblut_WH(nn) = 5/9 * rho_h2o * ds.tblut_retrieval.minTau * (ds.tblut_retrieval.minRe/1e6); % g/m^2


    % What is the true LWP
    lwp_inSitu(nn) = ds.GN_inputs.measurement.lwp;   % g/m^2
    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % ** Compute new updated LWP calc ***
    re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
        ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


    % define the z vector
    z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

    % define the z midpoint at each layer and normalize it!
    z_norm = z - z(1);
    z_norm_mid = (diff(z_norm)/2 + z_norm(1:end-1));


    % The radius input is defined as [r_start, r_end, r_step].
    % where r_step is the interval between radii values (used only for
    % vectors of radii). A 0 tells the code there is no step. Finally, the
    % radius values have to be in increasing order.
    ext_bulk_coeff_per_LWC = zeros(length(re_profile), 1);

    for rr = 1:length(re_profile)

        mie_radius = [re_profile(rr), re_profile(rr), 0];    % microns

        size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

        % Create a mie file
        [input_filename, output_filename] = write_mie_file('MIEV0', 'water',...
            mie_radius, 500, size_distribution, 'verbose', rr, round(re_profile(rr), 4), mie_folder_path);

        % run the mie file
        [~] = runMIE(mie_folder_path, input_filename,output_filename, which_computer);

        % Read the output of the mie file
        [mie,~,~] = readMIE(mie_folder_path, output_filename);

        ext_bulk_coeff_per_LWC(rr) = mie.Qext;       % km^-1 / (cm^3 / m^3)

    end


    % ** Assuming liquid water content increases linearly with depth **
    z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
    dz_km = z(2) - z(1);           % kilometers

    %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
    slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

    % solve for the linear liquid water content profile
    %lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
    lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


    lwp_newCalc(nn) = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    acpw_retrieval(nn) = ds.GN_outputs.retrieval(end,end);    % mm

    % What is the true LWP
    acpw_inSitu(nn) = ds.GN_inputs.measurement.actpw;   % mm
    % -------------------------------------------------------
    % -------------------------------------------------------


    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    tauC_retrieval(nn) = ds.GN_outputs.retrieval(3,end);    %

    % What is the true LWP
    tauC_inSitu(nn) = ds.GN_inputs.measurement.tau_c;   %
    % -------------------------------------------------------
    % -------------------------------------------------------



end

% -------------------------------------------------------
% Compute statistics!!

% Let's compute the root-mean-square percent error
rms_err_lwp_hyperspectral = 100 * sqrt( mean( (1 - lwp_retrieval./lwp_inSitu).^2 ));  % percent
rms_err_lwp_tblut = 100 * sqrt( mean( (1 - lwp_tblut./lwp_inSitu).^2 ));  % percent
rms_err_lwp_tblut_WH = 100 * sqrt( mean( (1 - lwp_tblut_WH./lwp_inSitu).^2 ));  % percent

% using the new LWP calc!
rms_err_lwp_hyperspectral_newCalc = 100 * sqrt( mean( (1 - lwp_newCalc./lwp_inSitu).^2 ));  % percent


% Let's compute the average percent difference
avg_percent_diff_newCacl = mean( abs( 100 .* (1 - lwp_newCalc./lwp_inSitu) ));
avg_percent_diff_tblut = mean( abs( 100 .* (1 - lwp_tblut./lwp_inSitu) ));
avg_percent_diff_tblutWH = mean( abs( 100 .* (1 - lwp_tblut_WH./lwp_inSitu) ));







%% Plot true color image from EMIT with overlapping footprints from the Aqua instruments


clear variables
% add libRadTran libraries to the matlab path
addLibRadTran_paths;
scriptPlotting_wht;


% Define EMIT Data locations and LibRadTran paths

folder_paths = define_EMIT_dataPath_and_saveFolders(2);
which_computer = folder_paths.which_computer;


% Would you like to print status updates and/or the libRadtran error file?


plot_figures = true;

save_figures = false;


% Define the folder of the coincident data set between EMIT and Aqau

% ---------------------------------------------
% ---------- PICK COINCIDENT DATA SET  --------
% ---------------------------------------------


% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    folder_paths.coincident_dataPath = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    folder_paths.coincident_dataFolder = '2024-09-12/';

elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    folder_paths.coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/'];

    % folder_paths.coincident_dataFolder = '2024-09-12/';

    % folder_paths.coincident_dataFolder = '2024_05_17-T1835/';

    % EMIT pixels masked out
    % folder_paths.coincident_dataFolder = '2023_9_16_T191106_2/';

    % 11 Pixels with H less than 1.6
    folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';

    % folder_paths.coincident_dataFolder = '2023_9_16_T191130_1/';

end








% Open Aqau data and look for overlapping pixels between EMIT and Aqua that meet certain criteria

criteria.cld_phase = 'water';
criteria.cld_cvr = 1;   % cloud fraction
criteria.cld_tau_min = 3;   % cloud optical depth
criteria.cld_tau_max = 30;   % cloud optical depth
criteria.H = 0.1;         % horizontal inhomogeneity index





% TODO: Add temporal information to compute time difference between pixels
% Will need:
% - EMIT pixel acquisition time (if available in emit.radiance structure)
% - MODIS pixel acquisition time (already available: modis.EV1km.pixel_time_UTC)
% - Then compute: overlap.time_difference_seconds(nn) = abs(emit_time(idx_emit) - modis_time(idx_modis))

[overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
    criteria, plot_figures);

% ** If there aren't any pixels found ... **
% Increase the horizontal inhomogeneity index

while isempty(overlap_pixels.modis.linear_idx) == true

    disp([newline, 'No overlaping pixels that meet defined criteria. Increasing H index....', newline])
    criteria.H = criteria.H + 0.25;         % horizontal inhomogeneity index

    % recompute
    [overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
        criteria, plot_figures);

    if isempty(overlap_pixels.modis.linear_idx) == false

        % print the H value used
        disp([newline, 'Horizontal Inhomogeneity Index - H = ', num2str(criteria.H), newline])

    end

end




% Plot all three swaths

if plot_figures == true

    figure; geoscatter(modis.geo.lat(:), modis.geo.long(:), 10, reshape(modis.cloud.effRadius17,[],1),'.');
    hold on; geoscatter(emit.radiance.geo.lat(:), emit.radiance.geo.long(:), 10, 'r.')
    hold on; geoscatter(airs.geo.Latitude(:), airs.geo.Longitude(:), 10, 'c.')
    hold on; geoscatter(amsr.geo.Latitude(:), amsr.geo.Longitude(:), 10, 'k.')

end




% Plot the pixel footprints on the Earth to see the overlap
% Add an RGB true color image for context

if plot_figures == true

    clear options
    % options.use_radiance = false;
    % options.rgb_image_type = 'modis';
    % [rgb_img, rgb_lat, rgb_lon] = create_modis_true_color(modis, options);

    options.convert_to_reflectance = false;
    options.rgb_image_type = 'emit';
    [rgb_img, rgb_lat, rgb_lon, band_indices] = create_emit_true_color(emit, options);


    options.show_rgb = true;
    options.rgb_image = rgb_img;
    options.rgb_lat = rgb_lat;
    options.rgb_lon = rgb_lon;



    % ** Plot with RGB Image **
    % fig = plot_instrument_footprints(modis, emit, amsr, overlap_pixels, options);
    % fig1 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    % [fig1, ax1] = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);
    % [fig1, ax1] = plot_instrument_footprints_4(modis, emit, airs, amsr, overlap_pixels, options);


    % Values for 2023_9_16_T191130_1 scene
    % options.latlim = [-25.7, -25.46];  % Only show -30° to -20° latitude
    % options.lonlim = [-71.15, -70.8];  % Only show -75° to -65° longitude

    % Values for 2023_9_16_T191118_1 scene
    options.latlim = [-25.8, -25.5];  % Only show -30° to -20° latitude
    options.lonlim = [-71.53, -71];  % Only show -75° to -65° longitude

    [fig1a, ax1a] = plot_instrument_footprints_4(modis, emit, [], amsr, overlap_pixels, options);


    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    % save .fig file
    if save_figures==true

        if strcmp(which_computer,'anbu8374')==true
            error(['Where do I save the figure?'])
        elseif strcmp(which_computer,'andrewbuggee')==true
            folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
        end
        f = gcf;
        saveas(f,[folderpath_figs,'EMIT Scene with MODIS context - ', folder_paths.coincident_dataFolder(1:end-1), '.fig']);


        % save .png with 500 DPI resolution
        % remove title
        ax1.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene with MODIS context - ', folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

        f = gcf;
        ax1a.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene with MODIS context - Zoom In version - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

    end
    % -------------------------------------
    % -------------------------------------


    % ** Plot without RGB Image **
    options.show_rgb = false;
    % fig2 = plot_instrument_footprints_2(modis, emit, amsr, overlap_pixels, options);
    fig2 = plot_instrument_footprints_3(modis, emit, airs, amsr, overlap_pixels, options);

    % ** Paper Worthy **
    % -------------------------------------
    % ---------- Save figure --------------
    if save_figures==true

        % save .fig file
        if strcmp(which_computer,'anbu8374')==true
            error(['Where do I save the figure?'])
        elseif strcmp(which_computer,'andrewbuggee')==true
            folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
        end
        % remove title
        ax1.Title.String = '';
        f = gcf;
        saveas(f,[folderpath_figs,'EMIT Scene and Aqua instrument overlap without MODIS context - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.fig']);


        % save .png with 400 DPI resolution
        % remove title
        ax1.Title.String = '';
        exportgraphics(f,[folderpath_figs,...
            'EMIT Scene and Aqau instrument overlap without MODIS context - ',...
            folder_paths.coincident_dataFolder(1:end-1), '.png'],'Resolution', 500);

    end
    % -------------------------------------
    % -------------------------------------

end



%%

% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);



%% Plot EMIT retrievals!

% clear variables





% Load simulated measurements
if strcmp(whatComputer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(whatComputer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % define the folder where the coincident data is stored
    data_path.stored_retrievals = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/Droplet_profile_retrievals/',...
        'Paper_2/Droplet_profile_retrievals_take1/'];

    addpath(data_path.stored_retrievals)



    % ------------------------------------------
    % **********  2023_9_16_T191130  ***********
    % ------------------------------------------
    % 2023_9_16_T191130 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_1_ran-on-08-Jan-2026_rev1.mat';

    % 2023_9_16_T191130 pixel 2
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191130_pixel_2_ran-on-09-Jan-2026_rev1.mat';



    % ------------------------------------------
    % **********  2023_9_16_T191118  ***********
    % ------------------------------------------
    % 2023_9_16_T191118 pixel 1
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_1_ran-on-08-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 2 -
    data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_2_ran-on-09-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 3 -
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_3_ran-on-09-Jan-2026_rev1.mat';

    % 2023_9_16_T191118 pixel 4 -
    % data_path.retrieval_name = '285bands_EMIT_dropRetrieval_2023_9_16_T191118_pixel_4_ran-on-09-Jan-2026_rev1.mat';


end

% make sure the MODIS, AIRS and AMSR-E data match the EMIT data used for
% the above retrieval
if strcmp(folder_paths.coincident_dataFolder(1:end-3),...
        extractBetween(data_path.retrieval_name, "Retrieval_", "_pixel")) == false

    error([newline, 'The MODIS, AIRS and AMSR-E data dont match the EMIT data time!', newline])

end

load([data_path.stored_retrievals, data_path.retrieval_name])

% what pixel number is this?
pixel_num_2Plot = str2double(extractBetween(data_path.retrieval_name, "pixel_", "_ran"));

% Make plot of the retrieved profile

% plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)
fig3 = plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR_perPixel(GN_outputs, GN_inputs, modis, [], amsr, pixel_num_2Plot);


% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
if strcmp(which_computer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(which_computer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig3,[folderpath_figs,'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1),...
    '_pixel_num_', num2str(pixel_num_2Plot),'.fig']);


% save .png with 500 DPI resolution
% remove title
exportgraphics(fig3,[folderpath_figs,...
    'EMIT Retrieval with MODIS and AMSR comparisons - ', folder_paths.coincident_dataFolder(1:end-1),...
    '_pixel_num_', num2str(pixel_num_2Plot), '.png'],'Resolution', 500);
% -------------------------------------
% -------------------------------------







%% Compute EMIT-Aqua retrieval statistics - take 1 - retrievals stored on CURC



clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_1/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';

    % % define the folder where retrievals are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_3/'];


    % define the folder where retrievals are located
    % *** 2/13/2026 - Retrieval with overlapping EMIT/Aqua data
    %          
    % folder_paths.retrieval = '/Users/andrewbuggee/MATLAB-Drive/EMIT/overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_5';


    % define the folder where retrievals are located
    % *** 2/15/2026 - Retrieval with overlapping EMIT/Aqua data ***
    % !! 37 retrieval !!
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
    %     'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7'];


    % define the folder where retrievals are located
    % *** 2/23/2026 - Retrieval with overlapping EMIT/Aqua data ***
    % !! 672 retrievals !!
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_12'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == false

        idx_2delete = [idx_2delete, nn];


    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];



con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3


% store the effective radii
re_modis = zeros(size(filenames_retrieval));
re_top = zeros(size(filenames_retrieval));
re_base = zeros(size(filenames_retrieval));

% store the LWP retrieval
% store the TBLUT LWP estimate
% store the TBLUT LWP estimate with Wood-Hartmann adjustement
% store the in-situ LWP measurement
lwp_retrieval = zeros(size(filenames_retrieval));
lwp_modis = zeros(size(filenames_retrieval));
lwp_modis_WH = zeros(size(filenames_retrieval));
lwp_newCalc = zeros(size(filenames_retrieval));

% store the LWP uncertainties for the MODIS retrieval
lwp_modis_err = zeros(size(filenames_retrieval));
lwp_modis_WH_err = zeros(size(filenames_retrieval));

% store the LWP estimate from AMSR measurements and their associated
% uncertainty
lwp_amsr_err = zeros(size(filenames_retrieval));
lwp_amsr = zeros(size(filenames_retrieval));

% store the ACPW retrieval
% store the true ACPW used in the forward model - measured by....
acpw_retrieval = zeros(size(filenames_retrieval));
acpw_modis = zeros(size(filenames_retrieval));
acpw_modis_corrected = zeros(size(filenames_retrieval));
acpw_airs = zeros(size(filenames_retrieval));

% store the AIRS acpw retrieval uncertainty estiamte
acpw_airs_err = zeros(size(filenames_retrieval));


% store the optical depth retrieval
% store the in-situ measured optical depth
tauC_retrieval = zeros(size(filenames_retrieval));
tauC_modis = zeros(size(filenames_retrieval));

% Store the MODIS optical thickness retrieval uncertainty
tauC_modis_err = zeros(size(filenames_retrieval));



idx_2delete_round2 = [];



for nn = 1:length(filenames_retrieval)

    clear ds emit modis airs

    ds = load([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name]);


    % Check to see if the retrieval converged
    if isfield(ds, 'GN_outputs')==false

        % skip this file
        idx_2delete_round2 = [idx_2delete_round2, nn];
        continue

    end


    % ----------------------------------------
    % *** Extract the EMIT pixel number ***
    % ----------------------------------------

    if strcmp(folder_paths.retrieval, ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7']) == true

        pixel_num = str2double(extractBetween([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name],...
            'pixel_', '_'));

    elseif strcmp(folder_paths.retrieval, ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_12']) == true

        pixel_num = 1;

    end




    
    % ----------------------------------------
    % *** Load MODIS, AIRS and AMSR-E data ***
    % ----------------------------------------

    % Load EMIT data
    % [emit, ~] = retrieveEMIT_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load Aqua/MODIS Data
    [modis, ~] = retrieveMODIS_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load AIRS data
    airs = readAIRS_L2_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load AMSR-E/2 data
    amsr = readAMSR_L2_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);
    % ----------------------------------------

    % ----------------------------------------
    % Remove data that is not needed
    % ----------------------------------------
    % emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

    modis = remove_unwanted_modis_data(modis, ds.overlap_pixels.modis);

    airs = remove_unwanted_airs_data(airs, ds.overlap_pixels.airs);

    amsr = remove_unwanted_amsr_data(amsr, ds.overlap_pixels.amsr);



    % Compute the above cloud precipitable water from AIRS data
    airs = convert_AIRS_prof_2_mass_density(airs, atm_data_directory,...
        pixel_num, ds.overlap_pixels, [], false, ds.GN_inputs.RT.z_topBottom(1)*1e3);



    % ** use the AIRS measurement closest to EMIT **
    unique_airs_pix = unique(ds.overlap_pixels.airs.linear_idx);
    unique_pix_idx_airs = zeros(1, length(ds.overlap_pixels.airs.linear_idx));
    for xx = 1:length(unique_pix_idx_airs)

        unique_pix_idx_airs(xx) = find(unique_airs_pix==ds.overlap_pixels.airs.linear_idx(xx));

    end


    % ** use the MODIS measurement closest to EMIT **
    unique_modis_pix = unique(ds.overlap_pixels.modis.linear_idx);
    unique_pix_idx_modis = zeros(1, length(ds.overlap_pixels.modis.linear_idx));
    for xx = 1:length(unique_pix_idx_modis)

        unique_pix_idx_modis(xx) = find(unique_modis_pix==ds.overlap_pixels.modis.linear_idx(xx));

    end



    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the MODIS retrieved effective radius
    re_modis(nn) = modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) );  % microns

    % store the radius at cloud top and base from the hyperspectral
    % retrieval
    re_top(nn) = ds.GN_outputs.retrieval(1,end);                     % microns
    re_base(nn) = ds.GN_outputs.retrieval(2,end);                    % microns


    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the liquid water paths
    lwp_retrieval(nn) = ds.GN_outputs.LWP;    % g/m^2

    % compute the LWP estimate using the TBLUT retrieval
    lwp_modis(nn) = (2 * rho_h2o *...
        (modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) )/1e6) *...
        modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) ) )/3; % g/m^2

    % compute the uncertainty of LWP - the uncertaities of effective radius
    % and optical depth add in quadrature and you must account for the
    % change in LWP with respect to each variable
    % 
    lwp_modis_err(nn) = 2/3 * rho_h2o * sqrt( (modis.cloud.optThickness17(unique_pix_idx_modis(pixel_num)) *...
        modis.cloud.effRadius17(unique_pix_idx_modis(pixel_num))/1e6 * 0.01*modis.cloud.effRad_uncert_17(unique_pix_idx_modis(pixel_num)) )^2 +...
        (modis.cloud.effRadius17(unique_pix_idx_modis(pixel_num))/1e6 *...
        modis.cloud.optThickness17(unique_pix_idx_modis(pixel_num)) * 0.01*modis.cloud.optThickness_uncert_17(unique_pix_idx_modis(pixel_num)) )^2 );

    % ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
    lwp_modis_WH(nn) = 5/9 * rho_h2o *...
        (modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) )/1e6) *...
        modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) ); % g/m^2

    % compute the uncertainty of LWP - the uncertaities of effective radius
    % and optical depth add in quadrature and you must account for the
    % change in LWP with respect to each variable
    % !! MODIS retrieval uncertainties are listed as percents !!
    lwp_modis_WH_err(nn) = 5/9 * rho_h2o * sqrt( (modis.cloud.optThickness17(unique_pix_idx_modis(pixel_num)) *...
        modis.cloud.effRadius17(unique_pix_idx_modis(pixel_num))/1e6 * 0.01*modis.cloud.effRad_uncert_17(unique_pix_idx_modis(pixel_num)) )^2 +...
        (modis.cloud.effRadius17(unique_pix_idx_modis(pixel_num))/1e6 *...
        modis.cloud.optThickness17(unique_pix_idx_modis(pixel_num)) * 0.01*modis.cloud.optThickness_uncert_17(unique_pix_idx_modis(pixel_num)) )^2 );

    % store the AMSR-E LWP estimate
    lwp_amsr(nn) = amsr.cloud.LiquidWaterPath;  % g/m^2
    % store the LWP error
    lwp_amsr_err(nn) = amsr.cloud.ErrorLWP;     % g/m^2

    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % ** Compute new updated LWP calc ***
    re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
        ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


    % define the z vector
    z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

    % define the z midpoint at each layer and normalize it!
    z_norm = z - z(1);
    z_norm_mid = (diff(z_norm)/2 + z_norm(1:end-1));


    % The radius input is defined as [r_start, r_end, r_step].
    % where r_step is the interval between radii values (used only for
    % vectors of radii). A 0 tells the code there is no step. Finally, the
    % radius values have to be in increasing order.
    ext_bulk_coeff_per_LWC = zeros(length(re_profile), 1);

    for rr = 1:length(re_profile)

        mie_radius = [re_profile(rr), re_profile(rr), 0];    % microns

        size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

        % Create a mie file
        [input_filename, output_filename] = write_mie_file('MIEV0', 'water',...
            mie_radius, 500, size_distribution, 'verbose', rr, round(re_profile(rr), 4), mie_folder_path);

        % run the mie file
        [~] = runMIE(mie_folder_path, input_filename,output_filename, which_computer);

        % Read the output of the mie file
        [mie,~,~] = readMIE(mie_folder_path, output_filename);

        ext_bulk_coeff_per_LWC(rr) = mie.Qext;       % km^-1 / (cm^3 / m^3)

    end


    % ** Assuming liquid water content increases linearly with depth **
    z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
    dz_km = z(2) - z(1);           % kilometers

    %slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
    slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

    % solve for the linear liquid water content profile
    %lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
    lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


    lwp_newCalc(nn) = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % store above cloud preciptiable water
    acpw_retrieval(nn) = ds.GN_outputs.retrieval(end,end);    % mm

    % What is the MODIS retrieved LWP
    acpw_modis(nn) = modis.vapor.col_nir( unique_pix_idx_modis(pixel_num) ) * 10;   % mm
    acpw_modis_corrected(nn) = modis.vapor.col_nir_corrected( unique_pix_idx_modis(pixel_num) ) * 10;   % mm

    % store the ACPW estimate from AIRS retrievals
    acpw_airs(nn) = airs.H2O.acpw_using_assumed_CTH( unique_pix_idx_airs(pixel_num) );  % mm

    % store the ACPW AIRS uncertainty estimate
    acpw_airs_err(nn) = airs.H2O.acpw_using_assumed_CTH_sigma;    % mm
    % -------------------------------------------------------
    % -------------------------------------------------------


    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the retrieved optical depth
    tauC_retrieval(nn) = ds.GN_outputs.retrieval(3,end);    %

    % What is the MODIS optical depth
    tauC_modis(nn) = modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) );

    % Store the optical thickness retrieval uncertainty
    tauC_modis_err(nn) = modis.cloud.optThickness_uncert_17(unique_pix_idx_modis(pixel_num));
    % -------------------------------------------------------
    % -------------------------------------------------------



end


% Remove indices with no converged solution
lwp_newCalc(idx_2delete_round2) = [];
lwp_modis(idx_2delete_round2) = [];
lwp_modis_WH(idx_2delete_round2) = [];
lwp_amsr(idx_2delete_round2) = [];

lwp_modis_err(idx_2delete_round2) = [];
lwp_modis_WH_err(idx_2delete_round2) = [];


acpw_retrieval(idx_2delete_round2) = [];
acpw_modis(idx_2delete_round2) = [];
acpw_modis_corrected(idx_2delete_round2) = [];
acpw_airs(idx_2delete_round2) = [];

acpw_airs_err(idx_2delete_round2) = [];

tauC_retrieval(idx_2delete_round2) = [];
tauC_modis(idx_2delete_round2) = [];
tauC_modis_err(idx_2delete_round2) = [];

re_modis(idx_2delete_round2) = [];
re_top(idx_2delete_round2) = [];

% -------------------------------------------------------
% Compute statistics!!

% % Let's compute the root-mean-square percent error
% rms_err_lwp_hyperspectral = 100 * sqrt( mean( (1 - lwp_retrieval./lwp_inSitu).^2 ));  % percent
% rms_err_lwp_tblut = 100 * sqrt( mean( (1 - lwp_tblut./lwp_inSitu).^2 ));  % percent
% rms_err_lwp_tblut_WH = 100 * sqrt( mean( (1 - lwp_tblut_WH./lwp_inSitu).^2 ));  % percent
% 
% % using the new LWP calc!
% rms_err_lwp_hyperspectral_newCalc = 100 * sqrt( mean( (1 - lwp_newCalc./lwp_inSitu).^2 ));  % percent


% Let's compute the average percent difference for LWP
avg_percent_LWP_diff_newCacl_MODIS = mean( abs( 100 .* (1 - lwp_newCalc./lwp_modis) ));
avg_percent_LWP_diff_newCalc_MODIS_WH = mean( abs( 100 .* (1 - lwp_newCalc./lwp_modis_WH) ));

avg_percent_LWP_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - lwp_newCalc./lwp_modis) ));
avg_percent_LWP_diff_newCalc_MODIS_WH_noAbs = mean( ( 100 .* (1 - lwp_newCalc./lwp_modis_WH) ));

% Compute the average percent difference for LWP between AMSR-E and the
% hyperspectral retrieval
idx_amsr_NOT_nan = find(~isnan(lwp_amsr));
avg_percent_LWP_diff_newCacl_AMSR_noAbs = mean( ( 100 .* (1 - lwp_newCalc(idx_amsr_NOT_nan)./...
    lwp_amsr(idx_amsr_NOT_nan)) ));


% Let's compute the average percent difference for ACPW
avg_percent_ACPW_diff_newCacl_MODIS = mean( abs( 100 .* (1 - acpw_retrieval./acpw_modis) ));
avg_percent_ACPW_diff_newCalc_AIRS = mean( abs( 100 .* (1 - acpw_retrieval./acpw_airs) ));

avg_percent_ACPW_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - acpw_retrieval./acpw_modis) ));
avg_percent_ACPW_diff_newCacl_MODIS_with_NIRcorrection_noAbs = mean( ( 100 .* (1 - acpw_retrieval./acpw_modis_corrected) ));
avg_percent_ACPW_diff_newCalc_AIRS_noAbs = mean( ( 100 .* (1 - acpw_retrieval./acpw_airs) ));

avg_percent_ACPW_diff_MODIS_AIRS_noAbs = mean( ( 100 .* (1 - acpw_airs./acpw_modis) ));


% Let's compute the average percent difference for optical thickness
avg_percent_tau_diff_newCacl_MODIS = mean( abs( 100 .* (1 - tauC_retrieval./tauC_modis) ));

avg_percent_tau_diff_newCacl_MODIS_noAbs = mean( ( 100 .* (1 - tauC_retrieval./tauC_modis) ));

% Compute the average retrieval bias between MODIS retrieved effective
% radius and the hyperspetral retrieval of radius at cloud top
re_rTop_avg_percent_diff = mean( ( 100 .* (1 - re_top./re_modis) ));


















%% Plot true and predicted EMIT spectrum with Mean relative error



clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/anbu8374/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_12/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    coincident_dataPath = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/Batch_Scripts/Paper-2/coincident_EMIT_Aqua_data/southEast_pacific/'];

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    atm_data_directory = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/';

    % define the folder where retrievals are located
    % *** 2/23/2026 - Retrieval with overlapping EMIT/Aqua data ***
    % !! 672 retrievals !!
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_12'];



elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end


% First step through all files in the directory and remove invisble files
% or directories

% Grab filenames in drive
filenames_retrieval = dir(folder_paths.retrieval);
idx_2delete = [];
for nn = 1:length(filenames_retrieval)

    if contains(filenames_retrieval(nn).name, "EMIT_dropRetrieval", "IgnoreCase", true) == false

        idx_2delete = [idx_2delete, nn];


    end

end

% delete rows that don't have retrieval filenames
filenames_retrieval(idx_2delete) = [];



con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3



% Define an optical depth range to plot
tauC_min = 5;
tauC_max = 6;



for nn = 1:length(filenames_retrieval)

    clear ds

    ds = load([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name]);


    % Check to see if the retrieval converged
    if isfield(ds, 'GN_outputs')==false

        % skip this file
        continue

    end

    % check that the optical depth is between the min and max values
    if (ds.GN_outputs.retrieval(3,end) < tauC_min) || (ds.GN_outputs.retrieval(3,end) > tauC_max)

        % skip this file
        continue

    end

   


    % ----------------------------------------
    % *** Extract the EMIT pixel number ***
    % ----------------------------------------

    if strcmp(folder_paths.retrieval, ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_7']) == true

        pixel_num = str2double(extractBetween([filenames_retrieval(nn).folder, '/', filenames_retrieval(nn).name],...
            'pixel_', '_'));

    elseif strcmp(folder_paths.retrieval, ['/Users/andrewbuggee/MATLAB-Drive/EMIT/',...
        'overlapping_with_Aqua/Droplet_profile_retrievals/Paper_2/take_12']) == true

        pixel_num = 1;

    end




    
    % ----------------------------------------
    % *** Load MODIS and EMIT data ***
    % ----------------------------------------

    % Load EMIT data
    [emit, ~] = retrieveEMIT_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);

    % Load Aqua/MODIS Data
    [modis, ~] = retrieveMODIS_data([coincident_dataPath, ds.folder_paths.coincident_dataFolder]);



    % ----------------------------------------

    % ----------------------------------------
    % Remove data that is not needed
    % ----------------------------------------
    emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

    modis = remove_unwanted_modis_data(modis, ds.overlap_pixels.modis);






    % ** use the MODIS measurement closest to EMIT **
    unique_modis_pix = unique(ds.overlap_pixels.modis.linear_idx);
    unique_pix_idx_modis = zeros(1, length(ds.overlap_pixels.modis.linear_idx));
    for xx = 1:length(unique_pix_idx_modis)

        unique_pix_idx_modis(xx) = find(unique_modis_pix==ds.overlap_pixels.modis.linear_idx(xx));

    end



    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the MODIS retrieved effective radius
    re_modis(nn) = modis.cloud.effRadius17( unique_pix_idx_modis(pixel_num) );  % microns

    % store the radius at cloud top and base from the hyperspectral
    % retrieval
    re_top(nn) = ds.GN_outputs.retrieval(1,end);                     % microns
    re_base(nn) = ds.GN_outputs.retrieval(2,end);                    % microns


    % -------------------------------------------------------
    % -------------------------------------------------------




    % -------------------------------------------------------
    % -------------------------------------------------------
    % store the retrieved optical depth
    tauC_retrieval(nn) = ds.GN_outputs.retrieval(3,end);    %

    % What is the MODIS optical depth
    tauC_modis(nn) = modis.cloud.optThickness17( unique_pix_idx_modis(pixel_num) );

    % Store the optical thickness retrieval uncertainty
    tauC_modis_err(nn) = modis.cloud.optThickness_uncert_17(unique_pix_idx_modis(pixel_num));
    % -------------------------------------------------------
    % -------------------------------------------------------








    % -------------------------------------------------------
    % -------------------------------------------------------
    % Grab the true measurement spectrum and the predicted spectrum
    
    % -------------------------------------------------------
    % -------------------------------------------------------



end
















%% Create plot showing hyperspectral retrieval of LWP with AMSR retrieval of LWP and its uncertainty

% Define the linewidth
ln_wdth = 2;

% define the font size
fnt_sz = 20;

% define marker width for the circle
circ_size = 30;

% define the equation font size
eq_fnt_sz = 16;


fig1 = figure;

% *** Plot retrieved LWP vs MODIS LWP ***
subplot(2,3,1)
errorbar(lwp_modis, lwp_newCalc,...
    [],[], lwp_modis_err, lwp_modis_err,...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_newCalc), 1.1 * max(lwp_newCalc)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(lwp_modis, lwp_newCalc, 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('MODIS LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_lwp = sprintf('y = %.3fx + %.3f', p_lwp(1), p_lwp(2));
text(0.05, 0.95, eq_str_lwp, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])
% plot legend
legend('', 'one-to-one', 'linear fit', 'Interpreter','latex',...
    'Position', [0.194507210511621 0.1673125 0.13906449167352 0.0699999999999998], 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')




% *** Plot hyperspectral LWP vs MODIS LWP w/ Wood-hartmann adjustment ***
subplot(2,3,2)
errorbar(lwp_modis_WH, lwp_newCalc,...
    [],[], lwp_modis_WH_err, lwp_modis_WH_err,...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_newCalc), 1.1 * max(lwp_newCalc)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(lwp_modis_WH, lwp_newCalc, 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('MODIS LWP_{WH} ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tblut = sprintf('y = %.3fx + %.3f', p_lwp(1), p_lwp(2));
text(0.05, 0.95, eq_str_tblut, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])





% *** Plot hyperspectral LWP vs AMSR LWP ***
subplot(2,3,3)
errorbar(lwp_amsr(idx_amsr_NOT_nan), lwp_newCalc(idx_amsr_NOT_nan),...
    [],[], lwp_amsr_err(idx_amsr_NOT_nan), lwp_amsr_err(idx_amsr_NOT_nan),...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_newCalc(idx_amsr_NOT_nan)), 1.1 * max(lwp_newCalc(idx_amsr_NOT_nan))];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(lwp_newCalc(idx_amsr_NOT_nan), lwp_amsr(idx_amsr_NOT_nan), 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('AMSR LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tblut_wh = sprintf('y = %.3fx + %.3f', p_lwp(1), p_lwp(2));
text(0.05, 0.95, eq_str_tblut_wh, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])





% *** Plot Tau ***
subplot(2,3,4)
errorbar(tauC_modis, tauC_retrieval,...
    [],[], tauC_modis_err, tauC_modis_err,...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(tauC_retrieval), 1.1 * max(tauC_retrieval)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_tau = polyfit(tauC_modis, tauC_retrieval, 1);
tau_fit = polyval(p_tau, ax_lim);
plot(ax_lim, tau_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('MODIS $\tau_c$', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral $\tau_c$', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tau = sprintf('y = %.3fx + %.3f', p_tau(1), p_tau(2));
text(0.05, 0.95, eq_str_tau, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])






% MODIS ACPW_err according to https://atmosphere-imager.gsfc.nasa.gov/products/water-vapor
acpw_modis_err = 0.1; 

% *** Plot MODIS ACPW versus hyperspectral retrieval ***
subplot(2,3,5)

errorbar(acpw_modis, acpw_retrieval,...
    [],[], acpw_modis*acpw_modis_err, acpw_modis*acpw_modis_err,...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(acpw_retrieval), 1.1 * max(acpw_retrieval)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_acpw = polyfit(acpw_modis, acpw_retrieval, 1);
acpw_fit = polyval(p_acpw, ax_lim);
plot(ax_lim, acpw_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('MODIS $IWV_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral $IWV_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_acpw = sprintf('y = %.3fx + %.3f', p_acpw(1), p_acpw(2));
text(0.05, 0.95, eq_str_acpw, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])



% *** Plot AIRS ACPW vs hyperspectral retrieval ***
subplot(2,3,6)

errorbar(acpw_airs, acpw_retrieval,...
    [],[], acpw_airs_err, acpw_airs_err,...
    '.','MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'),...
    'LineWidth', 1);
hold on
% plot a one to one line
ax_lim = [0.9 * min(acpw_retrieval), 1.1 * max(acpw_retrieval)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_acpw = polyfit(acpw_modis, acpw_retrieval, 1);
acpw_fit = polyval(p_acpw, ax_lim);
plot(ax_lim, acpw_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('AIRS $IWV_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Hyperspectral $IWV_{ac}$ ($mm$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_acpw = sprintf('y = %.3fx + %.3f', p_acpw(1), p_acpw(2));
text(0.05, 0.95, eq_str_acpw, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])





set(gcf,'Position',[0 0 1350 750])



% ** Paper Worthy **
% -------------------------------------
% ---------- Save figure --------------
% save .fig file
% if strcmp(whatComputer,'anbu8374')==true
%     error(['Where do I save the figure?'])
% elseif strcmp(whatComputer,'andrewbuggee')==true
%     folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
% end
% saveas(fig1,[folderpath_figs,'One-to-one comparison between retrieval of LWP,',...
%     'TauC and ACPW against the True values.fig']);
% 
% 
% % save .png with 500 DPI resolution
% % remove title
% exportgraphics(fig1,[folderpath_figs,'One-to-one comparison between retrieval of LWP,',...
%     'TauC and ACPW against the True values.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------



