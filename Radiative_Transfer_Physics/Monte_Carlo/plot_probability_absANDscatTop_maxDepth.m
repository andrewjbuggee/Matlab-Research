% Make a plot of either the probability or the PDF of the conditional
% probability: If photon n is scatterd out the top, the probability
% of reaching a maximum depth tau is...

% probability_str can either be
%   (1) 'probability' = which computes the number of counts in some bin and
%   divides it by the total number of photon
%   (2) 'pdf' = which computes the probability distribution function by
%   computing the ratio of the number of counts in a bin to the product of
%   the bin width and the total number of photons

function [] = plot_probability_absANDscatTop_maxDepth(inputs, final_state, photon_tracking, probability_str)


% ------------------------------------------------
% **************** Check Inputs ******************
% ------------------------------------------------

if strcmp(probability_str,'probability')==false && strcmp(probability_str,'pdf')==false

    error([newline, 'I dont recognize the histrogram normalization input',newline])
end



% ------------------------------------------------
% ------------------------------------------------


% First select those photons that were scattered out the top

index_scatter_out_top = final_state.scatter_out_top_INDEX;
num_photons_scatter_out_top = final_state.scatter_out_top;

[scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_edges] = ...
    histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);


% Next select those photons that were absorbed

index_absorbed = final_state.absorbed_INDEX;
num_photons_absorbed = final_state.absorbed;

[absorbed_maxDepth_PDF, absorbed_maxDepth_PDF_edges] = ...
    histcounts(photon_tracking.maxDepth(index_absorbed),'Normalization',probability_str);


% Make plot

if inputs.N_layers==1

    figure;
    plot(scatter_out_top_maxDepth_PDF,...
        scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2)
    hold on
    plot(absorbed_maxDepth_PDF,...
        absorbed_maxDepth_PDF_edges(1:end-1) + diff(absorbed_maxDepth_PDF_edges)/2)
    set(gca, 'YDir','reverse')
    grid on; grid minor
    xlabel('$P(\tau)$','Interpreter','latex');
    ylabel('$\tau$','Interpreter','latex')

    if strcmp(probability_str,'probability')==true
    title({'Conditional Probability of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')

    elseif strcmp(probability_str,'pdf')==true
        title({'Conditional PDF of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')
    end

    set(gcf, 'Position',[0 0 1000 630])

    dim = [0.685 0.5 0 0];

    texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
        ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
        ['$g$ = ', num2str(inputs.g)],...
        ['$r$ = ', num2str(inputs.layerBoundaries), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
        ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
    t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
    t.Color = 'black';
    t.FontSize = 25;
    t.FontWeight = 'bold';
    t.EdgeColor = 'black';
    t.FitBoxToText = 'on';

    legend(['Scattered out top, ', num2str(num_photons_scatter_out_top),' photons'],...
        ['Absorbed, ', num2str(num_photons_absorbed),' photons'],...
        'Interpreter','latex','Location','best');


else

    % There is more than one layer! For now lets leave out the ssa and
    % asymmetry parameter values

    figure;
    plot(scatter_out_top_maxDepth_PDF,...
        scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2)
    hold on
    plot(absorbed_maxDepth_PDF,...
        absorbed_maxDepth_PDF_edges(1:end-1) + diff(absorbed_maxDepth_PDF_edges)/2)
    set(gca, 'YDir','reverse')
    grid on; grid minor
    xlabel('$P(\tau)$','Interpreter','latex');
    ylabel('$\tau$','Interpreter','latex')

    if strcmp(probability_str,'probability')==true
    title({'Conditional Probability of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')

    elseif strcmp(probability_str,'pdf')==true
        title({'Conditional PDF of photons reaching a max depth of $\tau$'},...
        'Interpreter','latex')
    end

    set(gcf, 'Position',[0 0 1000 630])

    dim = [0.685 0.5 0 0];

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

    legend(['Scattered out top, ', num2str(num_photons_scatter_out_top),' photons'],...
        ['Absorbed, ', num2str(num_photons_absorbed),' photons'],...
        'Interpreter','latex','Location','best');


end