%% Plot five variables - truth vs retrieved with the one-to-one line

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------



elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';

    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];

    % % define the folder where retrievals are located
    % folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
    %     'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_with_reProf_uncert_3/'];


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

% Define the linewidth
ln_wdth = 2;

% define the font size
fnt_sz = 20;

% define marker width for the circle
circ_size = 30;

% define the equation font size
eq_fnt_sz = 16;

con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3



% store all retrieved values to plot
LWP_retrieved = zeros(length(filenames_retrieval), 1);
LWP_true = zeros(length(filenames_retrieval), 1);
lwp_tblut = zeros(size(filenames_retrieval));
lwp_tblut_WH = zeros(size(filenames_retrieval));

% ** try the new LWP calc!
lwp_newCalc = zeros(length(filenames_retrieval), 1);

tau_retrieved = zeros(length(filenames_retrieval), 1);
tau_true = zeros(length(filenames_retrieval), 1);

acpw_retrieved = zeros(length(filenames_retrieval), 1);
acpw_true = zeros(length(filenames_retrieval), 1);

idx_nan = [];

% how many retrievals failed or ran out of time?
no_retrieval = 0;


fig1 = figure;

for nn = 1:length(filenames_retrieval)


    % Load the data from the file
    ds = load([filenames_retrieval(nn).folder, '/',filenames_retrieval(nn).name]);




    % first, plot the simulated profile values as two lines
    if isfield(ds, 'GN_inputs')==true

        LWP_retrieved(nn) = ds.GN_outputs.LWP;
        LWP_true(nn) = ds.GN_inputs.measurement.lwp;

        % compute the LWP estimate using the TBLUT retrieval
        lwp_tblut(nn) = (2 * rho_h2o * (ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2
    
        % ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
        lwp_tblut_WH(nn) = 5/9 * rho_h2o * ds.tblut_retrieval.minTau * (ds.tblut_retrieval.minRe/1e6); % g/m^2



        tau_retrieved(nn) = ds.GN_outputs.retrieval(3,end);
        tau_true(nn) = ds.GN_inputs.measurement.tau_c;


        acpw_retrieved(nn) = ds.GN_outputs.retrieval(4,end);
        acpw_true(nn) = ds.GN_inputs.measurement.actpw;

        % *** new calculation of LWP ***
        re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
            ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);

        tau_vec = flipud(linspace(0, ds.GN_outputs.retrieval(3,end), length(re_profile)+1 )');
        
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




    else

        idx_nan  = [idx_nan, nn];

        no_retrieval = no_retrieval+1;

    end


    hold on


end

% remove nans

LWP_retrieved(idx_nan) = [];
LWP_true(idx_nan) = [];
lwp_newCalc(idx_nan) = [];

tau_retrieved(idx_nan) = [];
tau_true(idx_nan) = [];

acpw_retrieved(idx_nan) = [];
acpw_true(idx_nan) = [];



% *** Plot retrieved LWP vs in-situ LWP ***
subplot(1,5,1)
plot(LWP_true, lwp_newCalc, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_newCalc), 1.1 * max(lwp_newCalc)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(LWP_true, lwp_newCalc, 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('True LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
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




% *** Plot TBLUT LWP vs in-situ LWP ***
subplot(1,5,2)
plot(LWP_true, lwp_tblut, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_tblut), 1.1 * max(lwp_tblut)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(LWP_true, lwp_tblut, 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('True LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('TBLUT LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tblut = sprintf('y = %.3fx + %.3f', p_lwp(1), p_lwp(2));
text(0.05, 0.95, eq_str_tblut, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])





% *** Plot TBLUT-WH LWP vs in-situ LWP ***
subplot(1,5,3)
plot(LWP_true, lwp_tblut_WH, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(lwp_tblut_WH), 1.1 * max(lwp_tblut_WH)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_lwp = polyfit(LWP_true, lwp_tblut_WH, 1);
lwp_fit = polyval(p_lwp, ax_lim);
plot(ax_lim, lwp_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('True LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('$TBLUT_{WH}$ LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tblut_wh = sprintf('y = %.3fx + %.3f', p_lwp(1), p_lwp(2));
text(0.05, 0.95, eq_str_tblut_wh, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])





% *** Plot Tau ***
subplot(1,5,4)
plot(tau_true, tau_retrieved, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(tau_retrieved), 1.1 * max(tau_retrieved)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_tau = polyfit(tau_true, tau_retrieved, 1);
tau_fit = polyval(p_tau, ax_lim);
plot(ax_lim, tau_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('True $\tau_c$', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved $\tau_c$', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% Add textbox with equation
eq_str_tau = sprintf('y = %.3fx + %.3f', p_tau(1), p_tau(2));
text(0.05, 0.95, eq_str_tau, 'Units', 'normalized', 'FontSize', eq_fnt_sz,...
    'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'k',...
    'Position', [0.284442970682927 0.0479339227309894 0])






% *** Plot ACPW ***
subplot(1,5,5)

plot(acpw_true, acpw_retrieved, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(acpw_retrieved), 1.1 * max(acpw_retrieved)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
% Compute and plot linear fit
p_acpw = polyfit(acpw_true, acpw_retrieved, 1);
acpw_fit = polyval(p_acpw, ax_lim);
plot(ax_lim, acpw_fit, 'k--', 'LineWidth', 1)
grid on; grid minor
xlabel('True ACPW ($mm$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved ACPW ($mm$)', 'Interpreter','latex', 'FontSize', fnt_sz)
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
if strcmp(whatComputer,'anbu8374')==true
    error(['Where do I save the figure?'])
elseif strcmp(whatComputer,'andrewbuggee')==true
    folderpath_figs = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Presentations_and_Papers/paper_2/saved_figures/';
end
saveas(fig1,[folderpath_figs,'One-to-one comparison between retrieval of LWP,',...
    'TauC and ACPW against the True values.fig']);


% save .png with 400 DPI resolution
% remove title
exportgraphics(fig1,[folderpath_figs,'One-to-one comparison between retrieval of LWP,',...
    'TauC and ACPW against the True values.jpg'],'Resolution', 500);
% -------------------------------------
% -------------------------------------