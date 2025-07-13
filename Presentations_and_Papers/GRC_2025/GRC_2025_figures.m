%% Figures for GRC Poster!!

% By Andrew John Buggee

%% Plot droplet profile retrieval using 35 wavelengths with different forward model assumptions
% along with the retrieval using 66 wavelengths and the retrieval of CWV

clear variables


% Determine which computer you're using
which_computer = whatComputer();

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/',...
        'HySICS/Droplet_profile_retrievals/'];


elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % ---- Define where the retrievals are stored ---
    folder_paths.HySICS_retrievals = ['/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/',...
        'Hyperspectral_Cloud_Retrievals/HySICS/Droplet_profile_retrievals/'];




elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------

end


% define the mat files for each retreival to plot
filenames = {'dropletRetrieval_HySICS_35bands_withNoise_10mm-totalCWV_sim-ran-on-11-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_15mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_20mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_35bands_withNoise_25mm-totalCWV_sim-ran-on-12-Jul-2025_rev1.mat',...
    'dropletRetrieval_HySICS_66bands_withNoise_cwvRetrieval_sim-ran-on-12-Jul-2025_rev1.mat'};




% Step through each file

% define the colors for each curve plotted
C = mySavedColors(61:66, 'fixed');

figure;


for nn = 1:length(filenames)


    % Load a data set
    ds = load([folder_paths.HySICS_retrievals, filenames{nn}]);

    if nn==1

        % create a droplet profile from simulated measurement inputs
        re_sim = create_droplet_profile2([GN_inputs.RT.r_top, GN_inputs.RT.r_bot],...
            GN_outputs.tau_vector, 'optical_depth', GN_inputs.model.profile.type);

        tau_sim = linspace(0, GN_inputs.RT.tau_c, 100);

        % first, plot the simulated profile
        plot(re_sim, tau_sim, 'Color', C(1,:),'LineStyle','-', 'LineWidth',3)
        hold on

    end



    plot(GN_outputs.re_profile, GN_outputs.tau_vector, C(nn+1,:),'LineStyle',':', 'LineWidth',3)

    hold on


end


% flip the axes
set(gca, 'YDir','reverse')


% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.02,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.02,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, abs(GN_inputs.RT.z(end) - GN_inputs.RT.z(1))])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left

grid on; grid minor

ylim([-0.1, 1.2*GN_inputs.RT.tau_c])


set(gcf,'Position',[0 0 1200 630])

%%

