% Make a plot of either the probability or the PDF of the conditional
% probability: If photon n is scatterd out the top, the probability 
% of reaching a maximum depth tau is...

% probability_str can either be
%   (1) 'probability' = which computes the number of counts in some bin and
%   divides it by the total number of photon
%   (2) 'pdf' = which computes the probability distribution function by
%   computing the ratio of the number of counts in a bin to the product of
%   the bin width and the total number of photons

function [] = plot_probability_scatterOutTop_maxDepth(inputs, final_state, photon_tracking, probability_str)


% ------------------------------------------------
% **************** Check Inputs ******************
% ------------------------------------------------

if strcmp(probability_str,'probability')==false && strcmp(probability_str,'pdf')==false

    error([newline, 'I dont recognize the histrogram normalization input',newline])
end



% ------------------------------------------------
% ------------------------------------------------


% First select on those photons that were scattered out the top

index_scatter_out_top = final_state.scatter_out_top_index;
num_photons_scatter_out_top = final_state.scatter_out_top;

[scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_edges] = ...
    histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);

figure;
plot(scatter_out_top_maxDepth_PDF,...
    scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2)
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')
title({'Probability of a photon that scattered out the top', 'reaching a max depth of $\tau$'},...
    'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])

dim = [0.7 0.75 0 0];

if inputs.N_layers>1 && length(unique(inputs.layerRadii))>1

    texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
        ['N layers = ', num2str(inputs.N_layers)],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
        ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
        ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};

elseif inputs.N_layers==1


    texBox_str = {[num2str(num_photons_scatter_out_top),' photons'],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
        ['$g$ = ', num2str(inputs.g)],...
        ['$r$ = ', num2str(inputs.layerRadii), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_upper_limit)]};

end

t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'black';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'black';
t.FitBoxToText = 'on';


end