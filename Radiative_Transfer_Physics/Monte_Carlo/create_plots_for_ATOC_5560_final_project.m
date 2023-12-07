%% Create Plots for 5560 Final Project


% By Andrew John Buggee

%% Recreate S.Platnick (2000) figure 5a

clear variables

solar_zenith_angle = 0:10:70;

filenames = {'2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_03-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_10.mat',...
    '2D_MC_03-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_20.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_30.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_40.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_50.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_60.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_20_Tau0_8_SZA_70.mat'};


% Do you want to plot the probability or a pd?
probability_str = 'pdf';


% Define a set of colors based on the number of files
C = rand(length(filenames), 3);


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

else
    error('I dont recognize this computer user name')
end


if strcmp(pwd,saved_simulations)==false
    cd(saved_simulations)
end



% Start figure
figure;

for nn = 1:length(filenames)


    % Load a simulation
    load(filenames{nn})



    % First select those photons that were scattered out the top

    index_scatter_out_top = final_state.scatter_out_top_INDEX;

    [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_edges] = ...
        histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);


    % Plot the conditional probability
    plot(scatter_out_top_maxDepth_PDF,...
        scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2)
    hold on



    % Create legend string
    legend_str{nn} = ['$\mu_0 = ',num2str(round(cosd(solar_zenith_angle(nn)),2)),'$'];




end


% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Create title
if strcmp(probability_str,'probability')==true
    title({'Conditional Probability of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')

elseif strcmp(probability_str,'pdf')==true
    title({'Conditional PDF of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')
end


% Create textbox with simulation properties

% Textbox
dim = [0.685 0.5 0 0];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
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
legend(legend_str,'Interpreter','latex','Location','best','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])






%% Recreate S.Platnick (2000) figure 4

clear variables


wavelength = [660, 1650, 2155, 3700];

filenames = {'2D_MC_04-Dec-2022_Wavelength_660_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_1650_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_3700_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat'};

% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = [0.96,0.42,0.65;...          % Bubble gum pink
    0.49,0.92,0.04;...          % Neon green
    0.80,0.79,0.85;...          % A pale grey
    0.14,0.96,0.93];            % A bright electric blue


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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




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
        if nn~=1
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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end


    % horizontal line width
    horizontal_linewidth = 3;
    line_font_size = 16;

    % Plot line of constant average tau for 1.6 microns

    yline(tau_avg(2),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
        ['Depth of retrieved $r_e$ for ',num2str(wavelength(2)/1e3),' $\mu m$'], 'Interpreter','latex',...
        'FontSize',line_font_size,'LabelVerticalAlignment','bottom')

    % Plot line of constant average tau for 2.2 microns

    yline(tau_avg(3),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
        ['Depth of retrieved $r_e$ for ',num2str(wavelength(3)/1e3),' $\mu m$'], 'Interpreter','latex',...
        'FontSize',line_font_size,'LabelVerticalAlignment','top')

    % Plot line of constant average tau for 3.7 microns

    yline(tau_avg(4),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
        ['Depth of retrieved $r_e$ for ',num2str(wavelength(4)/1e3),' $\mu m$'], 'Interpreter','latex',...
        'FontSize',line_font_size,'LabelVerticalAlignment','top')







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
legend(legend_str,'Interpreter','latex','Location','best','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])


clear variables




%% Create Plot that shows re, ssa, and g as a function of cloud optical depth in 3 panels

clear variables

load('2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_1000_Tau0_4_SZA_0.mat');


% Define axes tick label font size
gca_ax_font_size = 14;

% Define x label font size
x_label_font_size = 14;

% Define legend font size
legend_font_size = 16;

% ----  plot re in the first subfigure -----
figure;

% Define a color for this curve
C1 = mySavedColors(1,'fixed');

s1 = subplot(1,3,1);
plot(inputs.dropletProfile.re, inputs.dropletProfile.tau_layer_mid_points,'Color',C1, 'LineStyle','-','LineWidth',4)
grid on; grid minor
set(gca,'YDir','reverse')
ylabel('$\tau$','Interpreter','latex', 'FontSize',17);
xlabel('$r_e$ $(\mu m)$','Interpreter','latex', 'FontSize', x_label_font_size);


% Set the yaxis tick labels
set(gca,'YTick', (0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
set(gca,'YTickLabel',(0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
set(gca,'FontSize',gca_ax_font_size)



% ----  plot asymmetry parameter in the second subfigure -----



% Define 2 colors for these curves
C2 = mySavedColors(8,'fixed');
C1 = mySavedColors(9,'fixed');

s2 = subplot(1,3,2);
plot(inputs.g, inputs.dropletProfile.tau_layer_mid_points,'Color',C1,'LineStyle','-','LineWidth',2)
hold on
plot(inputs.g_avg, inputs.dropletProfile.tau_layer_mid_points,'Color',C2,'LineStyle','-','LineWidth',3)
legend('$g(r)$','$<g(r_e)>$', 'Interpreter','latex', 'Location','best','Fontsize',legend_font_size)
grid on; grid minor
set(gca,'YDir','reverse')
xlabel('Asymmetry Parameter','Interpreter','latex', 'FontSize', x_label_font_size);

% Create the figure title
title('Cloud Optical Properties','Interpreter','latex')

% Set the yaxis tick labels
set(gca,'YTick', (0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
set(gca,'YTickLabel',[])
set(gca,'FontSize',gca_ax_font_size)



% ----  plot single scattering albedo in the third subfigure -----


s3 = subplot(1,3,3);
plot(inputs.ssa, inputs.dropletProfile.tau_layer_mid_points,'Color',C1,'LineStyle','-','LineWidth',2)
hold on
plot(inputs.ssa_avg, inputs.dropletProfile.tau_layer_mid_points,'Color',C2,'LineStyle','-','LineWidth',3)
legend('$\tilde{\omega}(r)$','$<\tilde{\omega}(r_e)>$', 'Interpreter','latex', 'location','best',...
    'Fontsize',legend_font_size)
grid on; grid minor
set(gca,'YDir','reverse')
xlabel('Single Scattering Albedo','Interpreter','latex', 'FontSize', x_label_font_size);


% Set the yaxis tick labels
set(gca,'YTick', (0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
set(gca,'YTickLabel',[])
set(gca,'FontSize',gca_ax_font_size)


% set size of the entire figure
set(gcf,"Position", [0 0 1000 630])



% Arrange the position of each subfigure

% Now shift the positions of the subfigures
% The vector is defined as [X_position, y_position, Height, Width]
s1.Position = s1.Position.*[1 1 1 1];

s2.Position = s2.Position.*[0.9 1 1 1];

s3.Position = s3.Position.*[0.9 1 1 1];



%% Create Plot showing the Weighting Function for smaller and smaller ssa vales


clear variables

% Define the set of wavelengths used for this analysis
wavelength = [578, 710, 915, 1145, 1390];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_710_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_915_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1145_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1390_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat'};


% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)

        % Clear the file variables because they take up a lot of memory
        clear F_norm photon_tracking final_state inputs

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
        if nn~=1
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
        %legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];
        if nn==1
            legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),6),6)];
        elseif nn==2
            legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),5),6)];
        elseif nn==3
            legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),4),6)];
        elseif nn==4
            legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),3),6)];
        elseif nn==5
            legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),2),6)];
        end



    end




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
legend(legend_str,'Interpreter','latex','Location','best','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])


% Clear the file variables because they take up a lot of memory
clear F_norm photon_tracking final_state inputs



%% Create a Plot showing the integrated retrieved radius using the same wavelengths and files as above



clear variables

% Define the set of wavelengths used for this analysis
wavelength = [578, 710, 915, 1145, 1390];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_710_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_915_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1145_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1390_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat'};


% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];

% Define the color of the radius curve plot define
radius_color = [0.827450980392157, 0.368627450980392, 0.376470588235294];           %  A pleasing salmon red


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

else
    error('I dont recognize this computer user name')
end


if strcmp(pwd,saved_simulations)==false
    cd(saved_simulations)
end


% Define plot line widths
% horizontal line width
horizontal_linewidth = 4;
line_font_size = 18;


% Start figure
figure;


% Plot the raw PDF's

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
    % create an re vector that is the same length as our weighting
    % function
    new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
    re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
    re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
    tau_avg(nn) = interp1(re, new_tau,re_avg);

    % -------------------------------------------------------------
    % -------------------------------------------------------------

    % Plot line of constant average tau

    if nn==length(wavelength)
        % Attach a label to our constant horizontal line

        yline(tau_avg(nn),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k',...
            'Label',['Depth of retrieved $r_e$ for $\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),2),6),' $\mu m$'],...
            'Interpreter','latex','FontSize',line_font_size,'LabelVerticalAlignment','bottom')

    else

        yline(tau_avg(nn),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color',C(nn,:),...
            'Interpreter','latex','FontSize',line_font_size,'LabelVerticalAlignment','bottom')

    end

    hold on

    % Plot the droplet radii curve
    if nn==length(wavelength)
        plot(re, new_tau,'Color',radius_color, 'LineWidth',4)
    end



    % Create legend string
    %legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];
    if nn==1
        legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),6),6)];
    elseif nn==2
        legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),5),6)];
    elseif nn==3
        legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),4),6)];
    elseif nn==4
        legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),3),6)];
    elseif nn==5
        legend_str{nn} = ['$\tilde{\omega}_0 = $',num2str(round(inputs.ssa_avg(1),2),6)];
    end




end







% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$r_e$ $(nm)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

% Define y axis limits
%ylim([5.5, 9])

% Create title
title({'Column averaged $r_e$ retrieval'},'Interpreter','latex')


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

% Add the droplet profile to the legend string
legend_str{end+1} = '$r_e$ Droplet profile';

% Create Legend
legend(legend_str,'Interpreter','latex','Location','best','FontSize',22)


set(gcf, 'Position',[0 0 1000 630])


% Clear the file variables because they take up a lot of memory
clear F_norm photon_tracking final_state inputs





%% How do the Weighting functions change with optical depth?




clear variables

% Define the set of wavelengths used for this analysis
Tau = [2,4,8,16,32];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_2_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_4_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_32_SZA_0.mat'};


% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));



% Define axes font size
gca_ax_font_size = 12;
y_label_font_size = 35;
x_label_font_size = 18;
font_size_title = 35;


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

else
    error('I dont recognize this computer user name')
end


if strcmp(pwd,saved_simulations)==false
    cd(saved_simulations)
end


% Store the subplot axes object
s = cell(1, length(Tau));


% Start figure
figure;

if smooth_curves==false

    % Plot the raw PDF's

    for nn = 1:length(filenames)

        % Clear the file variables because they take up a lot of memory
        clear F_norm photon_tracking final_state inputs

        % Load a simulation
        load(filenames{nn})



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % Plot the conditional probability
        s{nn} = subplot(1, length(Tau), nn);
        plot(scatter_out_top_maxDepth_PDF,...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on



        % Create Legend
        legend(['$\tau_0 = $',num2str((Tau(nn)))],'Interpreter','latex','Location','best',...
            'FontSize',22)




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)

        % Clear the file variables because they take up a lot of memory
        clear F_norm photon_tracking final_state inputs

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
        if nn~=1
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
        s{nn} = subplot(1,length(Tau),nn);
        plot(f(scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2),...
            scatter_out_top_maxDepth_PDF_tau_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_tau_edges)/2, 'Color',C(nn,:))
        hold on


        % Set the yaxis tick labels
        set(gca,'YTick', (0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
        set(gca,'YTickLabel',(0:inputs.tau_y_upper_limit/8:inputs.tau_y_upper_limit))
        set(gca,'FontSize',gca_ax_font_size)

        % Define y limits
        ylim([0, inputs.tau_y_upper_limit]);

        % Switch y axis direction
        set(gca, 'YDir','reverse')
        grid on; grid minor


        % set x axis label
        xlabel('$P(\tau)$','Interpreter','latex', 'FontSize',x_label_font_size);

        if nn==1
            % only make th y axis label on the first subfigure
            ylabel('$\tau$','Interpreter','latex', 'FontSize',y_label_font_size)


        end


        % Create title for the center subfigure only
        if nn==3
            title('Weighting Functions for different Optical Thicknesses', 'Interpreter','latex',...
                'FontSize',font_size_title)
        end


        % Arrange the position of each subfigure



        % Create Legend
        if nn==1 || nn==2
            legend(['$\tau_0 = $',num2str((Tau(nn)))],'Interpreter','latex','Location','southwest',...
                'FontSize',22)

        else
            legend(['$\tau_0 = $',num2str((Tau(nn)))],'Interpreter','latex','Location','southeast',...
                'FontSize',22)
        end





    end




end




% Now shift the positions of the subfigures
% The vector is defined as [X_position, y_position, Height, Width]
s{1}.Position = [0.0911923076923077 0.11 0.123739495798319 0.815];
s{2}.Position = [0.245866677440206 0.11 0.123739495798319 0.815];
s{3}.Position = [0.400541047188106 0.11 0.123739495798319 0.815];
s{4}.Position = [0.555215416936003 0.11 0.123739495798319 0.815];
s{5}.Position = [0.709889786683902 0.11 0.123739495798319 0.815];


% Plot the textbox withall relevant model parameters
% Create textbox with simulation properties

% Textbox
dim = [0.846538461538462 0.298055327884735 0.148878449660081 0.349563719734313];

texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';




% Set size of the entire figure
set(gcf, 'Position',[0 0 1300 630])


% Clear the file variables because they take up a lot of memory
clear F_norm photon_tracking final_state inputs




%% Create Bar Chart Plot with Final States for Multiple files with Different Tau values

clear variables

% Define the set of wavelengths used for this analysis
Tau = [2,4,8,16,32];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_2_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_4_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_32_SZA_0.mat'};


% Set up matrix
Y = zeros(length(Tau), 3);


% set up legend string
legend_str = cell(1, length(Tau));


% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];



for nn = 1:length(Tau)

    % Clear the file variables because they take up a lot of memory
    clear F_norm photon_tracking final_state inputs

    % Load a simulation
    load(filenames{nn},'final_state', 'inputs');

    Y(nn,:) = [final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]./inputs.N_photons;


    legend_str{nn} = ['$\tau_0 = $',num2str((Tau(nn)))];

end


figure;
X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});

b = bar(X,Y,'grouped');

% Set up bar labels
for nn = 1:length(Tau)
    bar_label_xLocation = b(nn).XEndPoints;
    bar_label_yLocation = b(nn).YEndPoints;
    bar_labels = string(round(b(nn).YData, 2));
    text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
    
    % set bar color
    b(nn).FaceColor = C(nn,:);
end



set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
set(gca, 'TitleFontSizeMultiplier',1.3)


legend(legend_str,'Interpreter','latex','Location','southeast',...
                'FontSize',22)


grid on; grid minor
ylabel('Percentage','Interpreter','latex')
title('Probability of each final state', 'Interpreter','latex')


set(gcf, 'Position',[0 0 1300 830])


ylim([0 1])


dim = [0.7 0.85 0 0];
texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';







%% Plot histogram of horizontal travel for multiple files that vary with tau

clear variables

% Define the set of wavelengths used for this analysis
Tau = [2,4,8,16,32];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_2_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_4_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_8_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_32_SZA_0.mat'};


% Set up matrix
Y = zeros(length(Tau), 3);


% set up legend string
legend_str = cell(1, length(Tau));



% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];


figure;

for nn = 1:length(Tau)

    % Clear the file variables because they take up a lot of memory
    clear F_norm photon_tracking final_state inputs

    % Load a simulation
    load(filenames{nn},'photon_tracking', 'inputs');
    
    s{nn} = subplot(1, length(Tau), nn);
    histogram(abs(photon_tracking.max_horizontal_position),'NumBins',100, 'EdgeColor','none',...
        'FaceColor',C(nn,:));
    hold on
    xlim([0, 150])

    legend(['$\tau_0 = $',num2str((Tau(nn)))],'Interpreter','latex','Location','northwest',...
                'FontSize',22)

    set(gca,'YScale', 'log')

    grid on; grid minor

    if nn == 1
        ylabel('Counts','Interpreter','latex')
    elseif nn==3
        xlabel('Horizontal Optical Depth','Interpreter','latex')
        title('Histogram of Maximum Horizontal Distance Travelled', 'Interpreter','latex')

    end


end


% Now shift the positions of the subfigures
% The vector is defined as [X_position, y_position, Height, Width]
s{1}.Position = [0.0911923076923077 0.11 0.123739495798319 0.815];
s{2}.Position = [0.245866677440206 0.11 0.123739495798319 0.815];
s{3}.Position = [0.400541047188106 0.11 0.123739495798319 0.815];
s{4}.Position = [0.555215416936003 0.11 0.123739495798319 0.815];
s{5}.Position = [0.709889786683902 0.11 0.123739495798319 0.815];





set(gcf, 'Position',[0 0 1300 730])


dim = [0.846538461538462 0.298055327884735 0.148878449660081 0.349563719734313];
texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';





%% Create my own version of the S.Platnick (2000) figure 4 weighting function plot

clear variables

% Define the set of wavelengths used for this analysis
%wavelength = [578, 915, 1145, 1390, 1450, 1550, 1750, 1850, 1950, 2155, 2250];

% filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_915_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1145_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1390_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1450_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1550_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1750_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1850_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_1950_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_04-Dec-2022_Wavelength_2155_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
%     '2D_MC_05-Dec-2022_Wavelength_2250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat'};


wavelength = [578, 1145, 1450, 1850, 1950, 2250];


filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1145_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1450_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1850_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_1950_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
    '2D_MC_05-Dec-2022_Wavelength_2250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat'};




% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = mySavedColors(1:length(wavelength), 'fixed');


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));

% Store horizontal labels
horizontal_label = cell(1, length(wavelength));


% keep track of horizontal line indexes
horizontal_line_index = [];


% Set a reasonable value for Number concentration

Nc = 125;           % total number of droplets per cubic centimeter


% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})


        % Compute the z component for the first wavelength
        if nn==1
            for ii = 1:length(inputs.dropletProfile.re)

                if ii==1
                    dz(ii) = inputs.dropletProfile.tau_layer_mid_points(ii)/...
                        (inputs.Qe_avg(ii) * pi * (inputs.layerRadii(ii) * 1e-6)^2 * Nc*(1e6));

                else
                    dz(ii) = (inputs.dropletProfile.tau_layer_mid_points(ii) - inputs.dropletProfile.tau_layer_mid_points(ii-1))/...
                        (inputs.Qe_avg(ii) * pi * (inputs.layerRadii(ii) * 1e-6)^2 * Nc*(1e6));

                end

            end

            total_z_depth = sum(dz);
        end



        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % -------------------------------------------------------------
        % Integrate the drolet profile with the weighting function to
        % get an average effective radius measured, and thus an average
        % optical depth.
        % -------------------------------------------------------------
        %if nn==1 || nn==7 || nn==8 || nn ==9 || nn==11
        if nn==1 || nn==5 || nn==6
            % create an re vector that is the same length as our weighting
            % function
            horizontal_line_index = [horizontal_line_index, nn];
            new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
            re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
            re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
            tau_avg(nn) = interp1(re, new_tau,re_avg);
            
            horizontal_label{nn} = ['Depth of retrieved $r_e$ for ',num2str(wavelength(nn)/1e3),' $\mu m$'];


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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end







end


% horizontal line width
horizontal_linewidth = 3;
line_font_size = 22;

for ii = 1:length(horizontal_line_index)


    if rem(ii,2)~=0

        % Plot line of constant average tau

        yline(tau_avg(horizontal_line_index(ii)),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
            horizontal_label{horizontal_line_index(ii)}, 'Interpreter','latex',...
            'FontSize',line_font_size,'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','right')


    elseif rem(ii,2)==0

        % Plot line of constant average tau

        yline(tau_avg(horizontal_line_index(ii)),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
            horizontal_label{horizontal_line_index(ii)}, 'Interpreter','latex',...
            'FontSize',line_font_size,'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','right')

    end


end






% Set axes tick label font size
set(gca,'FontSize',25)

% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex', 'FontSize',35)
ylabel('$\tau$','Interpreter','latex', 'FontSize',35)

% Create title
title({'Conditional probability of photons that scatter out cloud top',...
    'reaching a max depth of $\tau$'},'Interpreter','latex', 'FontSize',40)


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
legend(legend_str,'Interpreter','latex','Location','best','FontSize',25)


% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, total_z_depth])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30); 
yyaxis left



% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',...
    [0.04 0.1 0.0913076923076923 0.0422222222222221],...
    'String','Cloud Bottom',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.04 0.8 0.049769230769231 0.0666666666666665],...
    'String','Cloud Top',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');






set(gcf, 'Position',[0 0 1300 900])


clear variables






%% Create Platnick figure 4 with another set of wavelengths and a solar zenith angle of 45


clear variables


% Wavelengths used
% wavelength = [625, 875, 1050, 1250, 1400, 1590, 1625, 1850, 1900, 2150, 2250]; 



% filenames = {'2D_MC_08-Dec-2022_Wavelength_625_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_875_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1050_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1400_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1590_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1625_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1850_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_1900_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_2150_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
%     '2D_MC_08-Dec-2022_Wavelength_2250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat'};





% Wavelengths used
wavelength = [625, 1250, 1590, 1625, 1850, 1900, 2250]; 



filenames = {'2D_MC_08-Dec-2022_Wavelength_625_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_1250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_1590_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_1625_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_1850_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_1900_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat',...
    '2D_MC_08-Dec-2022_Wavelength_2250_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_45.mat'};




% Do you want to plot the probability of a set of PDF's?
probability_str = 'pdf';



% Do you want to smooth the raw PDF's?
smooth_curves = true;


% Define a set of colors based on the number of files
C = mySavedColors(1:length(wavelength), 'fixed');


% Store the number of photons from each simulation
legend_str = cell(1,length(filenames));

% Store horizontal labels
horizontal_label = cell(1, length(wavelength));


% keep track of horizontal line indexes
horizontal_line_index = [];

% Set a reasonable value for Number concentration

Nc = 125;           % total number of droplets per cubic centimeter

% Open folder where simulations are saved if it's not already open
% what computer are we using?


if strcmp(whatComputer,'anbu8374')

    saved_simulations = '/Users/anbu8374/Documents/MATLAB/Radiative_Transfer_Physics/Monte_Carlo/Monte_Carlo_Simulation_Results';



elseif strcmp(whatComputer,'andrewbuggee')

    saved_simulations = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Radiative_Transfer_Physics/',...
        'Monte_Carlo/Monte_Carlo_Simulation_Results'];

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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end



else

    % If this is true, we smooth each PDF to make a nice pretty plot, but
    % at the expense of loosing the PDF (the smoothed functions likely
    % won't integrate to 0)


    for nn = 1:length(filenames)


        % Load a simulation
        load(filenames{nn})



        % Compute the z component for the first wavelength
        if nn==1
            for ii = 1:length(inputs.dropletProfile.re)

                if ii==1
                    dz(ii) = inputs.dropletProfile.tau_layer_mid_points(ii)/...
                        (inputs.Qe_avg(ii) * pi * (inputs.layerRadii(ii) * 1e-6)^2 * Nc*(1e6));

                else
                    dz(ii) = (inputs.dropletProfile.tau_layer_mid_points(ii) - inputs.dropletProfile.tau_layer_mid_points(ii-1))/...
                        (inputs.Qe_avg(ii) * pi * (inputs.layerRadii(ii) * 1e-6)^2 * Nc*(1e6));

                end

            end

            total_z_depth = sum(dz);
        end

    


        % First select those photons that were scattered out the top

        index_scatter_out_top = final_state.scatter_out_top_INDEX;

        [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_tau_edges] = ...
            histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);



        % -------------------------------------------------------------
        % Integrate the drolet profile with the weighting function to
        % get an average effective radius measured, and thus an average
        % optical depth.
        % -------------------------------------------------------------
        if nn==1 || nn==6 || nn==7
            % create an re vector that is the same length as our weighting
            % function
            horizontal_line_index = [horizontal_line_index, nn];
            new_tau = linspace(inputs.dropletProfile.tau_layer_mid_points(1), inputs.dropletProfile.tau_layer_mid_points(end), length(scatter_out_top_maxDepth_PDF));
            re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, new_tau);
            re_avg = trapz(new_tau, re .* scatter_out_top_maxDepth_PDF);
            tau_avg(nn) = interp1(re, new_tau,re_avg);
            
            horizontal_label{nn} = ['Depth of retrieved $r_e$ for ',num2str(wavelength(nn)/1e3),' $\mu m$'];


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
        legend_str{nn} = ['$\lambda = ',num2str((wavelength(nn))),'$ nm'];




    end







end


% horizontal line width
horizontal_linewidth = 3;
line_font_size = 22;

for ii = 1:length(horizontal_line_index)


    if rem(ii,2)~=0

        % Plot line of constant average tau

        yline(tau_avg(horizontal_line_index(ii)),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
            horizontal_label{horizontal_line_index(ii)}, 'Interpreter','latex',...
            'FontSize',line_font_size,'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','right')


    elseif rem(ii,2)==0

        % Plot line of constant average tau

        yline(tau_avg(horizontal_line_index(ii)),'LineWidth',horizontal_linewidth, 'LineStyle',':','Color','k','Label',...
            horizontal_label{horizontal_line_index(ii)}, 'Interpreter','latex',...
            'FontSize',line_font_size,'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','right')

    end


end






% Set axes tick label font size
set(gca,'FontSize',25)

% Set up axes labels
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex', 'FontSize',35)
ylabel('$\tau$','Interpreter','latex', 'FontSize',35)

% Create title
title({'Conditional probability of photons that scatter out cloud top',...
    'reaching a max depth of $\tau$'},'Interpreter','latex', 'FontSize',40)


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
legend(legend_str,'Interpreter','latex','Location','best','FontSize',25)


% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, total_z_depth])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30); 
yyaxis left

% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',...
    [0.04 0.1 0.0913076923076923 0.0422222222222221],...
    'String','Cloud Bottom',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.04 0.8 0.049769230769231 0.0666666666666665],...
    'String','Cloud Top',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');




set(gcf, 'Position',[0 0 1300 900])


clear variables






%% Create Bar Chart Plot with Final States for Multiple files with Different Wavelenghts

clear variables

% Define the set of wavelengths used for this analysis
Wavelength = [578, 1390];

% Define the set of filenames to use
filenames = {'2D_MC_04-Dec-2022_Wavelength_578_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat',...
            '2D_MC_05-Dec-2022_Wavelength_1390_N-Photons_10000000_N-Layers_100_Tau0_16_SZA_0.mat'};

% Set up matrix
Y = zeros(length(Wavelength), 3);


% set up legend string
legend_str = cell(1, length(Wavelength));


% Define a set of colors based on the number of files
C = [0.979748378356085,   0.594896074008614,   0.117417650855806;...
    0.438869973126103,   0.262211747780845,   0.296675873218327;...
    0.111119223440599,   0.602843089382083,   0.318778301925882;...
    0.8492724501845559,  0.05437108503990062, 0.9681090252965144;...
    0.3563645953575846,  0.4380836048512262,  0.5147715889386915];



for nn = 1:length(Wavelength)

    % Clear the file variables because they take up a lot of memory
    clear F_norm photon_tracking final_state inputs

    % Load a simulation
    load(filenames{nn},'final_state', 'inputs');

    Y(nn,:) = [final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]./inputs.N_photons;


    legend_str{nn} = ['$\lambda_0 = $',num2str((Wavelength(nn))), ' $nm$'];

end


figure;
X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});

b = bar(X,Y,'grouped');

% Set up bar labels
for nn = 1:length(Wavelength)
    bar_label_xLocation = b(nn).XEndPoints;
    bar_label_yLocation = b(nn).YEndPoints;
    bar_labels = string(round(b(nn).YData, 2));
    text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
    
    % set bar color
    b(nn).FaceColor = C(nn,:);
end



set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
set(gca, 'TitleFontSizeMultiplier',1.3)


legend(legend_str,'Interpreter','latex','Location','southeast',...
                'FontSize',22)


grid on; grid minor
ylabel('Percentage','Interpreter','latex')
title('Probability of each final state', 'Interpreter','latex')


set(gcf, 'Position',[0 0 1300 630])


ylim([0 1])


dim = [0.7 0.85 0 0];
texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ['N layers = ', num2str(inputs.N_layers)],...
    ['$\tau_0$ = ',num2str(inputs.tau_y_upper_limit)],...
    ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
    ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
    ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';










