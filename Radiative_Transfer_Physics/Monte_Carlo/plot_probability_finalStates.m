% Plot a bar chart showing the proportions of photons ending in each final
% state

function [] = plot_probability_finalStates(final_state,inputs)


if inputs.N_layers==1

    % There is only one layer in our medium!


    figure;
    X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});
    b = bar(X,[final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]);
    set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
    set(gca, 'TitleFontSizeMultiplier',1.3)
    bar_label_xLocation = b.XEndPoints;
    bar_label_yLocation = b.YEndPoints;
    bar_labels = string(b.YData./inputs.N_photons);
    text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
    grid on; grid minor
    ylabel('Counts','Interpreter','latex')
    title('Probability of each final state', 'Interpreter','latex')
    set(gcf, 'Position',[0 0 1000 630])
    ylim([0 inputs.N_photons])


    dim = [0.7 0.85 0 0];
    texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
        ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
        ['$g$ = ', num2str(inputs.g)],...
        ['$r$ = ', num2str(inputs.layerRadii), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
        ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
    t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
    t.Color = 'black';
    t.FontSize = 25;
    t.FontWeight = 'bold';
    t.EdgeColor = 'black';
    t.FitBoxToText = 'on';



else

    % There is more than 1 layer in the modeled medium!
    % Let's ignore plotting the ssa and the asymmetry parameter for now

    figure;
    X = categorical({'Scatter out top','Scatter out Bottom','Absorbed'});
    b = bar(X,[final_state.scatter_out_top, final_state.scatter_out_bottom, final_state.absorbed]);
    set(gca,'TickLabelInterpreter', 'latex','FontSize',25);
    set(gca, 'TitleFontSizeMultiplier',1.3)
    bar_label_xLocation = b.XEndPoints;
    bar_label_yLocation = b.YEndPoints;
    bar_labels = string(b.YData./inputs.N_photons);
    text(bar_label_xLocation,bar_label_yLocation,bar_labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontSize',20,'FontWeight','bold','Interpreter','latex')
    grid on; grid minor
    ylabel('Counts','Interpreter','latex')
    title('Probability of each final state', 'Interpreter','latex')
    set(gcf, 'Position',[0 0 1000 630])
    ylim([0 inputs.N_photons])


    dim = [0.7 0.85 0 0];
    texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
        ['N layers = ', num2str(inputs.N_layers)],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
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


end