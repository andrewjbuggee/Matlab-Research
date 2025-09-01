%% Analysis of retrievals from simulated HySICS measurements with negligible radiometric uncertainty


% By Andrew John Buggee

%% Plot a curve of r-top versus r-bot for each wavelength

clear variables


% Access specific file or folder
filenames = dir('/Users/andrewbuggee/MATLAB-Drive/vza_7_vaz_210_sza_10_saz_91/*.mat');

% what are the free parameters?
r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);

num_files = num_rTop * num_rBot * num_tauC * num_tcpw;

% Step through each file

% For some simulated above cloud total column precipitable water amount,
% let's show the retrieval of r-top and r-bot for different optical depths.
% But somehow I want to show a comparison between the simulated values and
% the retrieved values for each independent variable.

% But that isn't how my data is arranged. I have a retrieval set where
% every profile has a radius at clod top of 10 microns. So I will have to
% step through the other variables for now.

% How do I look at the difference between the retrieved values and the true
% value for 4 different variables? Could I have panels for different values
% of the above cloud precipitable water?

% The y axis is the percent difference between the true r-top value and the
% retrieved

% The x axis is the percent difference between the true r-bot value and the
% retrieved

% Each curve represents a different optical depth

% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tcpw-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);



figure;


for nn = 1:length(filenames)

    if nn==1

        % make a 1-to-1 line
        plot(tau_c, tau_c, 'k-', 'LineWidth',1)

        grid on; grid minor
        hold on

    end


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        % acpw_true = [acpw_true, ds.GN_inputs.measurement.actpw];

        [~, idx_actpw] = min(abs(ds.GN_inputs.measurement.actpw - acpw_true));

        plot(ds.GN_inputs.measurement.tau_c, ds.GN_outputs.retrieval(3,end), '.', 'MarkerSize', 20,...
            'Color', C(idx_actpw, :))

        hold on

    end


end

legend(["one to one", "$acpw = $"+string((round(acpw_true, 2)))], 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

title('$r_{top} = 10$ simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)
ylabel('$\tau_c$ retrieved', 'Interpreter','latex', 'FontSize',30)
xlabel('$\tau_c$ true', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 750])


%% Plot 11 panels for each acpw amount and show the percent difference between the true and retrieved radius at cloud top and bottom


clear variables


% Access specific file or folder
filenames = dir('/Users/andrewbuggee/MATLAB-Drive/vza_7_vaz_210_sza_10_saz_91/*.mat');

% what are the free parameters?
r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);

num_files = num_rTop * num_rBot * num_tauC * num_tcpw;

% Step through each file

% For some simulated above cloud total column precipitable water amount,
% let's show the retrieval of r-top and r-bot for different optical depths.
% But somehow I want to show a comparison between the simulated values and
% the retrieved values for each independent variable.

% But that isn't how my data is arranged. I have a retrieval set where
% every profile has a radius at clod top of 10 microns. So I will have to
% step through the other variables for now.

% How do I look at the difference between the retrieved values and the true
% value for 4 different variables? Could I have panels for different values
% of the above cloud precipitable water?

% The y axis is the percent difference between the true r-top value and the
% retrieved

% The x axis is the percent difference between the true r-bot value and the
% retrieved

% Each curve represents a different optical depth

% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tauC-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);



figure;
for nn = 1:num_tcpw

    subplot(3,4,nn)
    grid on; grid minor
    hold on
    title(['$acpw = $', num2str(acpw_true(nn))], 'Interpreter','latex', 'FontSize',20)


    if nn==9

        ylabel('$r_{top}^{true}$ - $r_{top}^{retrieved}$', 'Interpreter','latex', 'FontSize',20)
        xlabel('$r_{bot}^{true}$ - $r_{bot}^{retrieved}$', 'Interpreter','latex', 'FontSize',20)

    end


end


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        idx_tau_c = ds.GN_inputs.measurement.tau_c==tau_c;

        [~, idx_actpw] = min(abs(ds.GN_inputs.measurement.actpw - acpw_true));

        subplot(3,4,idx_actpw)

        hold on

        diff_rTop = ds.GN_inputs.measurement.r_top - ds.GN_outputs.retrieval(1,end);  % microns
        diff_rBot = ds.GN_inputs.measurement.r_bot - ds.GN_outputs.retrieval(2,end);  % microns

        plot(diff_rBot, diff_rTop, '.', 'MarkerSize', 20,...
            'Color', C(idx_tau_c, :))

        hold on

    end


end

for nn = 1:num_tcpw

    subplot(3,4,nn)
    grid on; grid minor
    hold on
    title(['$acpw =$ ', num2str(acpw_true(nn)), ' $mm$'], 'Interpreter','latex', 'FontSize',20)


    if nn==9

        ylabel('$r_{top}^{true}$ - $r_{top}^{retrieved}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',25)
        xlabel('$r_{bot}^{true}$ - $r_{bot}^{retrieved}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',25)

    end


end

legend("$\tau_c = $"+string((round(tau_c, 2))), 'Position',...
    [0.78153780173431 0.0609496949905001 0.057275146484375 0.257175660160735],...
    'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

%title('$r_{top} = 10$ simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)


set(gcf,'Position',[0 0 1500 875])




%% Let's look at the true acpw and the retrieved acpw as a function of cloud optical depth




clear variables


% Access specific file or folder
filenames = dir('/Users/andrewbuggee/MATLAB-Drive/vza_7_vaz_210_sza_10_saz_91/*.mat');

% what are the free parameters?
r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);

num_files = num_rTop * num_rBot * num_tauC * num_tcpw;

% Step through each file

% For some simulated above cloud total column precipitable water amount,
% let's show the retrieval of r-top and r-bot for different optical depths.
% But somehow I want to show a comparison between the simulated values and
% the retrieved values for each independent variable.

% But that isn't how my data is arranged. I have a retrieval set where
% every profile has a radius at clod top of 10 microns. So I will have to
% step through the other variables for now.

% How do I look at the difference between the retrieved values and the true
% value for 4 different variables? Could I have panels for different values
% of the above cloud precipitable water?

% The y axis is the percent difference between the true r-top value and the
% retrieved

% The x axis is the percent difference between the true r-bot value and the
% retrieved

% Each curve represents a different optical depth

% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tauC-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);



figure;

% first plot the multispectral resutls
acpw_tblut = zeros(num_tcpw, num_tauC);
% Then plot the hyperspectral retrieval
acpw_hyper = zeros(num_tcpw, num_tauC);



for nn = 1:length(filenames)



    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        idx_tau_c = find(ds.GN_inputs.measurement.tau_c==tau_c);

        [~, idx_actpw] = min(abs(ds.GN_inputs.measurement.actpw - acpw_true));

        acpw_tblut(idx_actpw, idx_tau_c) = ds.acpw_retrieval.min_interpolated;

        acpw_hyper(idx_actpw, idx_tau_c) = ds.GN_outputs.retrieval(end,end);    % mm




    end


end


% plot the TBLUT estimate
subplot(1,2,1)

% plot a one-to-one line
plot(linspace(0,max(acpw_true)+3, 100), linspace(0,max(acpw_true)+3, 100), 'k', 'LineWidth', 1)
hold on
grid on; grid minor

for tt = 1:num_tauC

    plot(acpw_tblut(:, tt), acpw_true, '.-', 'Color', C(tt,:), 'Linewidth', 1, ...
        'MarkerSize', 35)
    hold on

end




title('$r_{top} = 10$ simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)
ylabel('$acpw$ true', 'Interpreter','latex', 'FontSize',30)
xlabel('$acpw$ TBLUT', 'Interpreter','latex', 'FontSize',30)

legend(["one-to-one", "$\tau_c = $"+string((round(tau_c, 2)))],...
    'Interpreter','latex', 'Location','northwest', 'FontSize', 20)


% plot the hyperspectral measurement
subplot(1,2,2)

% plot a one-to-one line
plot(linspace(0,max(acpw_true)+3, 100), linspace(0,max(acpw_true)+3, 100), 'k', 'LineWidth', 1)
hold on
grid on; grid minor

for tt = 1:num_tauC

    plot(acpw_hyper(:, tt), acpw_true, '.-', 'Color', C(tt,:), 'Linewidth', 1, ...
        'MarkerSize', 35)
    hold on

end




ylabel('$acpw$ true', 'Interpreter','latex', 'FontSize',30)
xlabel('$acpw$ hyperspectral', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 950 950])




%% Let's compare the true value to the retrieved for all four variables in individual panels

clear variables


% Access specific file or folder
filenames = dir('/Users/andrewbuggee/MATLAB-Drive/vza_7_vaz_210_sza_10_saz_91/*.mat');

% what are the free parameters?
r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];


% set empty arrays
rTop_true = [];
rTop_retrieved = [];

rBot_true = [];
rBot_retrieved = [];

tau_c_true = [];
tau_c_retrieved = [];

acpw_true = [];
acpw_retrieved = [];



for nn = 1:length(filenames)



    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        rTop_true = [rTop_true, ds.GN_inputs.measurement.r_top];
        rTop_retrieved = [rTop_retrieved, ds.GN_outputs.retrieval(1,end)];

        rBot_true = [rBot_true, ds.GN_inputs.measurement.r_bot];
        rBot_retrieved = [rBot_retrieved, ds.GN_outputs.retrieval(2,end)];

        tau_c_true = [tau_c_true, ds.GN_inputs.measurement.tau_c];
        tau_c_retrieved = [tau_c_retrieved, ds.GN_outputs.retrieval(3,end)];

        acpw_true = [acpw_true, ds.GN_inputs.measurement.actpw];
        acpw_retrieved = [acpw_retrieved, ds.GN_outputs.retrieval(4,end)];




    end


end




% % plot a one-to-one line
% plot(linspace(0,max(rTop_true)+3, 100), linspace(0,max(rTop_true)+3, 100), 'k', 'LineWidth', 1)
% hold on
% grid on; grid minor
% plot(rTop_retrieved, rTop_true, '.', 'Color', mySavedColors(64, 'fixed'),'MarkerSize', 35)
%
% ylabel('$r_{top}$ true', 'Interpreter','latex', 'FontSize',30)
% xlabel('$r_{top}$ retrieved', 'Interpreter','latex', 'FontSize',30)



figure;
% plot just three variables for now

subplot(1,3,1)

% plot a one-to-one line
plot(linspace(0,max(rBot_true)+4, 100), linspace(0,max(rBot_true)+4, 100), 'k', 'LineWidth', 1)
hold on
grid on; grid minor
plot(rBot_retrieved, rBot_true, '.', 'Color', mySavedColors(62, 'fixed'),'MarkerSize', 35)
ylim([0, max(rBot_true)+4])

% compute the linaer fit
% fitresult1 = fit(rBot_retrieved', rBot_true', 'linear');

ylabel('$r_{bot}$ true', 'Interpreter','latex', 'FontSize',30)
xlabel('$r_{bot}$ retrieved', 'Interpreter','latex', 'FontSize',30)



subplot(1,3,2)

% plot a one-to-one line
plot(linspace(0,max(tau_c_true)+4, 100), linspace(0,max(tau_c_true)+4, 100), 'k', 'LineWidth', 1)
hold on
grid on; grid minor
plot(tau_c_retrieved, tau_c_true, '.', 'Color', mySavedColors(63, 'fixed'),'MarkerSize', 35)
ylim([0, max(tau_c_true)+4])

ylabel('$\tau_{c}$ true', 'Interpreter','latex', 'FontSize',30)
xlabel('$\tau_{c}$ retrieved', 'Interpreter','latex', 'FontSize',30)




subplot(1,3,3)

% plot a one-to-one line
plot(linspace(0,max(acpw_true)+4, 100), linspace(0,max(acpw_true)+4, 100), 'k', 'LineWidth', 1)
hold on
grid on; grid minor
plot(acpw_retrieved, acpw_true, '.', 'Color', mySavedColors(64, 'fixed'),'MarkerSize', 35)
ylim([0, max(acpw_true)+4])


ylabel('$acpw$ true', 'Interpreter','latex', 'FontSize',30)
xlabel('$acpw$ retrieved', 'Interpreter','latex', 'FontSize',30)

set(gcf,'Position',[0 0 1250 950])


%% Let's show the absolute difference between the true and retrieved r-top as a function of optical depth
% Each panel will be for a fixed radius at cloud bottom, and markers will
% be colored by their assumed tcpw



clear variables


% Access specific file or folder
filenames = dir('/Users/andrewbuggee/MATLAB-Drive/vza_7_vaz_210_sza_10_saz_91/*.mat');

% what are the free parameters?
r_top = 10;
r_bot = 3:10;
tau_c = 5:3:29;
tcpw = 5:3:35;

acpw_true = [2.8932, 4.6291, 6.3650, 8.1010,...
    9.8369, 11.5728, 13.3088, 15.0447, 16.7806,...
    18.5165, 20.2525];

% length of each independent variable
num_rTop = length(r_top);
num_rBot = length(r_bot);
num_tauC = length(tau_c);
num_tcpw = length(tcpw);


% Step through each file


% define the colors for each curve plotted
C = mySavedColors(61:(61+num_tcpw-1), 'fixed');

% lgnd_str = cell(1, length(filenames) + 2);



figure;
for nn = 1:num_rBot

    subplot(2,4,nn)
    grid on; grid minor
    hold on
    title(['$r_{bot} = $', num2str(r_bot(nn)) ' $\mu m$'], 'Interpreter','latex', 'FontSize',20)


end


for nn = 1:length(filenames)


    % Load a data set
    ds = load([filenames(nn).folder, '/', filenames(nn).name]);

    if isfield(ds, 'GN_outputs')

        idx_rBot = find(ds.GN_inputs.measurement.r_bot==r_bot);

        tau_c_true = ds.GN_inputs.measurement.tau_c;

        [~, idx_actpw] = min(abs(ds.GN_inputs.measurement.actpw - acpw_true));

        subplot(2,4,idx_rBot)

        hold on

        diff_rTop = ds.GN_inputs.measurement.r_top - ds.GN_outputs.retrieval(1,end);  % microns

        plot(tau_c_true, diff_rTop, '.', 'MarkerSize', 20,...
            'Color', C(idx_actpw, :))

        hold on

    end


end

for nn = 1:num_rBot

    subplot(2,4,nn)
    hold on
    title(['$r_{bot} = $', num2str(r_bot(nn)) ' $\mu m$'], 'Interpreter','latex', 'FontSize',20)

    if nn==1

        ylabel('$r_{top}^{true}$ - $r_{top}^{retrieved}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)

    end

    if nn==5

        ylabel('$r_{top}^{true}$ - $r_{top}^{retrieved}$ ($\mu m$)', 'Interpreter','latex', 'FontSize',30)
        xlabel('$\tau_c$', 'Interpreter','latex', 'FontSize',30)

    end


    if nn>=6 || nn<=8

        xlabel('$\tau_c$', 'Interpreter','latex', 'FontSize',30)

    end


end

legend("$acpw = $"+string((round(acpw_true, 2))), 'Position',...
    [0.78153780173431 0.0609496949905001 0.057275146484375 0.257175660160735],...
    'Interpreter','latex', 'Location','northwest', 'FontSize', 20)

%title('$r_{top} = 10$ simulated HySICS measurements', 'Interpreter','latex', 'FontSize',30)


set(gcf,'Position',[0 0 1500 875])



%% I wan't to determine if there is a relationship between the percent different of some retrieved variable and the other three retrieved variables

% For each retrieved variable, I could compute the percent difference
% between the true and retrieved value. Then I could plot this three times
% against the three indpendent variables. There will be redundant values
% though, because for 