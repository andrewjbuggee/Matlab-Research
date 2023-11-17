% Make a plot of either the probability or the PDF of the max depth each
% photon reaches within our simulation

% probability_str can either be
%   (1) 'probability' = which computes the number of counts in some bin and
%   divides it by the total number of photon
%   (2) 'pdf' = which computes the probability distribution function by
%   computing the ratio of the number of counts in a bin to the product of
%   the bin width and the total number of photons

function [] = plot_probability_maxDepth(inputs, photon_tracking, probability_str)


% ------------------------------------------------
% **************** Check Inputs ******************
% ------------------------------------------------

if strcmp(probability_str,'probability')==false && strcmp(probability_str,'pdf')==false

    error([newline, 'I dont recognize the histrogram normalization input',newline])
end



% ------------------------------------------------
% ------------------------------------------------


% create a pdf of the probability a photon reaches a maximum penetration
% depth at some tau

[photon_tracking.maxDepth_PDF, photon_tracking.maxDepth_PDF_edges] = ...
    histcounts(photon_tracking.maxDepth,'Normalization',probability_str);

figure;
plot(photon_tracking.maxDepth_PDF,...
    photon_tracking.maxDepth_PDF_edges(1:end-1) + diff(photon_tracking.maxDepth_PDF_edges)/2)
set(gca, 'YDir','reverse')
grid on; grid minor
xlabel('$P(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')
title('Probability of a photon reaching a max depth of $\tau$', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1000 630])

dim = [0.7 0.85 0 0];

texBox_str = {[num2str(inputs.N_photons),' photons'],...
    ['$\lambda$ = ',num2str(inputs.wavelength), ' $nm$'],...
    ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
    ['$g$ = ', num2str(inputs.g)],...
    ['$r$ = ', num2str(inputs.radius), ' $\mu m$'],...
    ['$\tau_0$ = ', num2str(inputs.tau_upper_limit)],...
    ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
t.Color = 'white';
t.FontSize = 25;
t.FontWeight = 'bold';
t.EdgeColor = 'white';
t.FitBoxToText = 'on';


end