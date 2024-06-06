%% Create figures for IRS 2024 Presentation

% By Andrew John Buggee



%% Reflectance spectra example over cloudy pixel

clear variables

% --- Load EMIT Data ---

% -------------------------------------
% ------- PICK EMIT DATA SET  --------
% -------------------------------------

emitDataFolder = '17_Jan_2024_coast/';

% -------------------------------------


% Determine which computer you're using

% Find the folder where the mie calculations are stored
% find the folder where the water cloud files are stored.
if strcmp(whatComputer,'anbu8374')==true

    % ------ Folders on my Mac Desktop --------

    % Define the EMIT data folder path

    emitDataPath = '/Users/anbu8374/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP files will be saved
    folder2save.libRadTran_INP_OUT = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/'];

    % Define the folder path where the mat files of reflectances will be
    % saved
    folder2save.reflectance_calcs = emitDataPath;


elseif strcmp(whatComputer,'andrewbuggee')==true

    % ------ Folders on my Macbook --------

    % Define the EMIT data folder path

    emitDataPath = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/EMIT/EMIT_data/';


    % Define the folder path where all .INP and .OUT files will be saved
    folder2save.libRadTran_INP_OUT = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
        'Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/'];

    % Define the folder path where the mat files of reflectances will be
    % saved
    folder2save.reflectance_calcs = emitDataPath;




end

[emit,L1B_fileName] = retrieveEMIT_data([emitDataPath, emitDataFolder]);


% --- Define the pixels to use for the retrieval ---

% 17_Jan_2024_coast - optical depth of 6.6
pixels2use.row = 932;
pixels2use.col = 960;

% Grab the pixel indices
pixels2use = grab_pixel_indices(pixels2use, size(emit.radiance.measurements));



% --- Remove data that is not needed ---

emit = remove_unwanted_emit_data(emit, pixels2use);

inputs = create_emit_inputs_TBLUT(emitDataFolder, folder2save, emit);

% --- Define the spectral response function of EMIT for the desired Bands
% ---

emit.spec_response = create_EMIT_specResponse(emit, inputs);

% --- Define the solar source file name and read in the solar source data
% ---

% ********* IMPORTANT *************
% The source flux is integrated with the EMIT spectral response function

% define the source file using the input resolution
inputs = define_source_for_EMIT(inputs, emit);


% --- Convert radiance measurements to TOA reflectance for the desired
% pixels ---

emit = convert_EMIT_radiance_2_reflectance(emit, inputs);

% --- Create plot ---

figure;
plot(emit.radiance.wavelength, emit.reflectance.value, '.-', 'MarkerSize', 20,...
    'LineWidth',1, 'Color', mySavedColors(3, 'fixed'))
xlabel('Wavelength ($nm$)', Interpreter='latex', FontSize=30)
ylabel('Reflectance ($1/sr$)', Interpreter='latex', FontSize=30)
grid on; grid minor

% set figure size
set(gcf, 'Position', [0 0 1200 650])


%%         RAW WEIGHTING FUNCTIONS



% Plot weighting functions of the first 7 MODIS spectral channels
% These weighting functions were created using the VOCALS-REx data set from
% Nov-9-2023. The droplet profile and optical depth were modeled after the
% vertical profile sampled at 1.734 hrs after the plane took off. The SZA
% was set as the value measured by MODIS for the median pixel, the pixel
% found closest the C130 aircraft in the middle of it's ascent through the
% cloud.

% read contents of Monte Carlo Simulation folder
folder_contents = dir('/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results/');
filenames = {};

for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>5

        if strcmp(folder_contents(nn).name(1:17), '2D_MC_06-Jun-2024')==true

            filenames = [filenames, {folder_contents(nn).name}];

        end

    end

end

% plot as a pdf
probability_str = 'pdf';

% define the wavelengths as the changing variables
changing_variable = [499.771423339844
551.866699218750
611.462219238281
671.097534179688
768.055725097656
872.517578125000
1014.29498291016
1036.67773437500
1044.13830566406
1245.51538085938
1252.97241210938
1260.42834472656
1267.88330078125
1275.33923339844
1558.43286132813
1565.87658691406
1573.31933593750
1580.76208496094
1588.20495605469
1595.64672851563
1603.08862304688
1610.52954101563
1617.97045898438
1625.41040039063
1632.85131835938
1640.29028320313
1647.73034667969
1655.16943359375
1662.60742187500
2063.69653320313
2130.41748046875
2226.71948242188
2234.12329101563
2241.52685546875
2248.92968750000
2256.33276367188
2263.73461914063];

plot_probability_scat_top_maxDepth_with_changing_variable(filenames, probability_str ,changing_variable)


%%    SMOOTH WEIGHTING FUNCTIONS


clear variables


% read contents of Monte Carlo Simulation folder
folder_contents = dir('/Users/anbu8374/Documents/MATLAB/Matlab-Research/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results/');
filenames = {};

for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>5

        if strcmp(folder_contents(nn).name(1:17), '2D_MC_06-Jun-2024')==true

            filenames = [filenames, {folder_contents(nn).name}];


        end

    end

end


% ---------------------------------------------------------------------------------------


% Do you want to smooth the raw PDF's?
smooth_curves = true;

% Do you want to plot a horizontal line for the effective retrieved radius?
plot_retrieved_radius = false;

% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';

% Define a set of colors based on the number of files
C = mySavedColors(1:length(filenames), 'fixed');


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = ['/Users/anbu8374/Documents/MATLAB/Matlab-Research/',...
                    'Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results'];



elseif strcmp(whatComputer,'andrewbuggee')

    error([newline, 'Where is the new folder?', newline])

else
    error('I dont recognize this computer user name')
end


if strcmp(pwd,saved_simulations)==false
    cd(saved_simulations)
end



% Start figure
figure;

if smooth_curves==false

    % Plot the raw PDF's

    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % Plot the conditional probability
        plot(scatter_out_top_maxDepth_PDF,...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on



        % Create legend string
        legend_str{nn} = ['$\lambda = ',num2str(round(inputs.wavelength)),'$ nm'];




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % -------------------------------------------------------------
        % Integrate the drolet profile with the weighting function to
        % get an average effective radius measured, and thus an average
        % optical depth.
        % -------------------------------------------------------------
        if nn~=0
            % create an re vector that is the same length as our weighting
            % function
            new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
            re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
            re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
            tau_avg(nn) = interp1(re, new_tau,re_avg);


        end
        % -------------------------------------------------------------
        % -------------------------------------------------------------

        % Create smooth spline function
        f=fit((scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2)',scatter_out_top_maxDepth_PDF', 'smoothingspline','SmoothingParam',0.95);

        % Plot the conditional probability
        plot(f(scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2),...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on



        % Create legend string
        legend_str{nn} = ['$\lambda = ',num2str(round(inputs.wavelength)),'$ nm'];




    end

    if plot_retrieved_radius==true
        % horizontal line width
        horizontal_linewidth = 4;
        line_font_size = 23;

        for nn = 1:length(filenames)
            % Plot line of constant tau associated with retrieval depth
            yline(tau_avg(nn),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color',C(nn,:),'Label',...
                ['Depth of retrieved $r_e$ for ',num2str(round(inputs.wavelength)/1e3),' $\mu m$'], 'Interpreter','latex',...
                'FontSize',line_font_size,'LabelVerticalAlignment','middle')

        end

    end

%     % Plot line of constant average tau for 0.66 microns
% 
%     yline(tau_avg(1),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
%         ['Depth of retrieved $r_e$ for ',num2str(wavelength(1)/1e3),' $\mu m$'], 'Interpreter','latex',...
%         'FontSize',line_font_size,'LabelVerticalAlignment','bottom')
% 
%     % Plot line of constant average tau for 1.6 microns
% 
%     yline(tau_avg(2),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
%         ['Depth of retrieved $r_e$ for ',num2str(wavelength(2)/1e3),' $\mu m$'], 'Interpreter','latex',...
%         'FontSize',line_font_size,'LabelVerticalAlignment','bottom')
% 
%     % Plot line of constant average tau for 2.2 microns
% 
%     yline(tau_avg(3),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
%         ['Depth of retrieved $r_e$ for ',num2str(wavelength(3)/1e3),' $\mu m$'], 'Interpreter','latex',...
%         'FontSize',line_font_size,'LabelVerticalAlignment','top')
% 
%     % Plot line of constant average tau for 3.7 microns
% 
%     yline(tau_avg(4),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
%         ['Depth of retrieved $r_e$ for ',num2str(wavelength(4)/1e3),' $\mu m$'], 'Interpreter','latex',...
%         'FontSize',line_font_size,'LabelVerticalAlignment','top')







end



% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
title({'Conditional probability of photons that scatter out cloud top',...
    'reaching a max depth of $\tau$'},'Interpreter','latex')


% Create textbox with simulation properties

% Textbox
dim = [0.685 0.5 0 0];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';


% Create Legend
legend(legend_str,'Interpreter','latex','Location','northwest','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])


%clear variables

