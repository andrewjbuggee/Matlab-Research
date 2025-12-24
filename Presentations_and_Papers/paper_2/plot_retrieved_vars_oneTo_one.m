%% Plot three variables - truth vs retrieved with the one-to-one line

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


    % define the folder where the Vocals-Rex in-situ derived measurements
    % are located
    folder_paths.retrieval = ['/Users/andrewbuggee/MATLAB-Drive/HySICS/Droplet_profile_retrievals/',...
        'paper2_variableSweep/full_retrieval_logSpace_newCov_VR_meas_allBands_3/'];



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



% store all retrieved values to plot
LWP_retrieved = zeros(length(filenames_retrieval), 1);
LWP_true = zeros(length(filenames_retrieval), 1);

tau_retrieved = zeros(length(filenames_retrieval), 1);
tau_true = zeros(length(filenames_retrieval), 1);

acpw_retrieved = zeros(length(filenames_retrieval), 1);
acpw_true = zeros(length(filenames_retrieval), 1);

idx_nan = [];

% how many retrievals failed or ran out of time?
no_retrieval = 0;


figure;

for nn = 1:length(filenames_retrieval)


    % Load the data from the file
    ds = load([filenames_retrieval(nn).folder, '/',filenames_retrieval(nn).name]);

    


    % first, plot the simulated profile values as two lines
    if isfield(ds, 'GN_inputs')==true

        LWP_retrieved(nn) = ds.GN_outputs.LWP;
        LWP_true(nn) = ds.GN_inputs.measurement.lwp;



        tau_retrieved(nn) = ds.GN_outputs.retrieval(3,end);
        tau_true(nn) = ds.GN_inputs.measurement.tau_c;


        acpw_retrieved(nn) = ds.GN_outputs.retrieval(4,end);
        acpw_true(nn) = ds.GN_inputs.measurement.actpw;




    else

        idx_nan  = [idx_nan, nn];

        no_retrieval = no_retrieval+1;

    end


    hold on


end

% remove nans

LWP_retrieved(idx_nan) = [];
LWP_true(idx_nan) = [];

tau_retrieved(idx_nan) = [];
tau_true(idx_nan) = [];

acpw_retrieved(idx_nan) = [];
acpw_true(idx_nan) = [];



% *** Plot LWP ***
subplot(1,3,1)
plot(LWP_true, LWP_retrieved, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
ax_lim = [0.9 * min(LWP_retrieved), 1.1 * max(LWP_retrieved)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
grid on; grid minor
xlabel('True LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved LWP ($g/m^{2}$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% plot legend
legend('', 'one-to-one', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')


% *** Plot Tau ***
subplot(1,3,2)


plot(tau_true, tau_retrieved, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));

hold on
% plot a one to one line
subplot(1,3,2)
ax_lim = [0.9 * min(tau_retrieved), 1.1 * max(tau_retrieved)];
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
grid on; grid minor
xlabel('True $\tau_c$', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved $\tau_c$', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% plot legend
legend('', 'one-to-one', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')


% *** Plot ACPW ***
subplot(1,3,3)

plot(acpw_true, acpw_retrieved, '.',...
    'MarkerSize', circ_size, 'Color', mySavedColors(63, 'fixed'));
hold on
% plot a one to one line
subplot(1,3,3)
ax_lim = [0.9 * min(acpw_retrieved), 1.1 * max(acpw_retrieved)];
grid on; grid minor
plot(ax_lim, ax_lim, 'k-', 'LineWidth', 1)
xlabel('True ACPW ($mm$)', 'Interpreter','latex', 'FontSize',fnt_sz)
ylabel('Retrieved ACPW ($mm$)', 'Interpreter','latex', 'FontSize', fnt_sz)
xlim(ax_lim)
ylim(ax_lim)
% plot legend
legend('', 'one-to-one', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')



set(gcf,'Position',[0 0 950 750])