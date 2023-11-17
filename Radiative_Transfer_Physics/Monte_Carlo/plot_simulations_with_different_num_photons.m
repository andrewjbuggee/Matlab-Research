%% Create Plot to Compare F_down and F_up for the same simulation but with different number of photons

% Andrew Buggee
%%

function plot_simulations_with_different_num_photons(filenames)

if iscell(filenames)~=true
    error([newline,'Enter more than one filename using a celll structure',newline])
end


%% Loop through file names


% Define colors for F_up and F_down

% Make the usual plot but without the 2 stream curves

C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];


% Store the number of photons from each simulation
legend_str = cell(1,2*length(filenames));

% Cycle through different linestyles
linestyles_for_plot = {':','-','--','-.'};

% Start figure
figure;

for nn = 1:length(filenames)


    % Load a simulation
    load(filenames{nn})


    
    plot(F_norm.up,F_norm.binEdges(1:end-1)+diff(F_norm.binEdges)/2,'Color',C1,...
        'LineStyle',linestyles_for_plot{nn},'Linewidth',2)
    hold on
    plot(F_norm.down,F_norm.binEdges(1:end-1)+diff(F_norm.binEdges)/2, 'Color',C2,...
        'LineStyle',linestyles_for_plot{nn},'LineWidth',2)
    
    
    legend_str{2*nn-1} = ['$F_{\uparrow}/F_0$, $N_{photons} = 10^{',num2str(log10(inputs.N_photons)),'}$'];
    legend_str{2*nn} = ['$F_{\downarrow}/F_0$, $N_{photons} = 10^{',num2str(log10(inputs.N_photons)),'}$'];




end


% Set up axes labels
grid on; grid minor
set(gca,'YDir','reverse')
xlabel('$F/F_0$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex');
title({'Comparing Monte Carlo simulations with different', 'total photon counts'}, 'Interpreter','latex')




% Create textbox with simulation properties

% Textbox
dim = [0.7 0.85 0 0];

texBox_str = {['N layers = ', num2str(inputs.N_layers)],...
    ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
    ['$\mu_0$ = ',num2str(cosd(inputs.solar_zenith_angle))],...
    ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
    ['$g$ = ', num2str(inputs.g)],...
    ['$r$ = ', num2str(inputs.mie.radius(1)), ' $\mu m$'],...
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

