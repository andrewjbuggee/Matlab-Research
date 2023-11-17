%% Plot Probability of Scattering out the top after reaching a max depth of Tau for multiple files


% By Andrew John Buggee



%% Create Plot to Compare F_down and F_up for the same simulation but with different number of photons

% Andrew Buggee
%%

function plot_probability_scat_top_maxDepth_with_changing_variable(filenames, probability_str ,changing_variable)

if iscell(filenames)~=true
    error([newline,'Enter more than one filename using a celll structure',newline])
end


%% Loop through file names



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



    % ----------------------------- Create legend string --------------------------
    % for changing solar zenith angle
    %legend_str{nn} = ['$\mu_0 = ',num2str(round(cosd(changing_variable(nn)),2)),'$'];
    
    % for changing wavelength
    legend_str{nn} = ['$\lambda = $',num2str(round(changing_variable(nn))),' $\mu m$'];




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

texBox_str = {['$N_{photons} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
    ...['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\theta_0$ = ',num2str((inputs.solar_zenith_angle)), ' $^{\circ}$'],...
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

