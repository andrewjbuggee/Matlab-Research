%% Plot the EMIT retrieved droplet profile and optical depth for a single pixel

% INPUTS:



% By Andrew John Buggee

%%

function plot_EMIT_retrieved_vertProf(emit, retrieval, tblut_retrieval)


% determine the number of pixels
num_pixels = size(emit.reflectance.value, 2);

% Plot the Gauss-Newton Retrieval
for pp = 1:num_pixels

    figure;

    plot(retrieval.re_profile(:,pp), retrieval.tau_vector(:,pp), 'Color',...
        mySavedColors(1,'fixed'),'LineStyle',':', 'LineWidth',3)

    % flip y-axis and provide axes labels
    set(gca,'YDir','reverse')
    ylabel('$\tau$','interpreter','latex','FontSize',35);
    xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
    grid on; grid minor; hold on;

    % Plot the retrieval uncertainty of the radius at cloud top
    errorbar(retrieval.re_profile(1,pp), retrieval.tau_vector(1,pp), sqrt(retrieval.posterior_cov(1,1,pp)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the radius at cloud bottom
    errorbar(retrieval.re_profile(end,pp), retrieval.tau_vector(end,pp), sqrt(retrieval.posterior_cov(2,2,pp)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the optical depth
    errorbar(retrieval.re_profile(end,pp), retrieval.tau_vector(end,pp), sqrt(retrieval.posterior_cov(3,3,pp)),...
        'vertical', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)



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



    % Plot the emit TBLUT droplet estimate as a constant vertical line

    xl0 = xline(tblut_retrieval.minRe(pp),':',...
        ['$$r_{e} = $$',num2str(round(tblut_retrieval.minRe(pp), 1)), '$$\mu m$$'], 'Fontsize',22,...
        'Interpreter','latex','LineWidth',3,'Color', mySavedColors(9, 'fixed'));
    xl0.LabelVerticalAlignment = 'bottom';

    % Plot the emit optical depth TBLUT retrieval as a constant horizontal line
    yl0 = yline(tblut_retrieval.minTau(pp),':',...
        ['$$\tau_{c} = $$',num2str(round(tblut_retrieval.minTau(pp), 1))], 'Fontsize',22,...
        'Interpreter','latex','LineWidth',3,'Color', mySavedColors(9, 'fixed'));
    yl0.LabelVerticalAlignment = 'top';
    yl0.LabelHorizontalAlignment = 'left';


    % compute the LWP estimate using the TBLUT retrieval
    rho_liquid_water = 10^6;        % g/m^3

    lwp_emit_tblut = (2*rho_liquid_water*(tblut_retrieval.minRe(pp)/1e6) * tblut_retrieval.minTau(pp))/3; % g/m^2

    % grab the hypersepctral retrieval estimate of LWP
    retrieved_LWP = retrieval.LWP;        % g/m^2

    % Print this information on the figure

    dim = [.137 .35 .3 .3];
    str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_emit_tblut,1)),' $g/m^{2}$', newline,...
        '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];

    annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'FontWeight','bold');
    set(gcf,'Position',[0 0 1200 630])


    % Create a Legend with only the two black curves
    legend({'Hyperspectral retrieval', '', '', '', 'TBLUT $r_e$', 'TBLUT $\tau_c$'},...
        'Interpreter','latex', 'Location','northwest', 'FontSize', 25)

    % set figure size
    set(gcf,'Position',[0 0 1200 630])



end



end