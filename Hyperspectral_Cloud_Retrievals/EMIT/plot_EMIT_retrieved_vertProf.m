%% Plot the EMIT retrieved droplet profile and optical depth for a single pixel

% INPUTS:



% By Andrew John Buggee

%%

function plot_EMIT_retrieved_vertProf(GN_outputs, tblut_retrieval, GN_inputs)


C = mySavedColors(1:2, 'fixed');

% Plot the Gauss-Newton Retrieval

figure;

title('Retrieved droplet profile using EMIT', 'Interpreter','latex',...
    'FontSize', 26)

plot(GN_outputs.re_profile, GN_outputs.tau_vector', 'Color',...
    C(1,:),'LineStyle',':', 'LineWidth',3)

% flip y-axis and provide axes labels
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
grid on; grid minor; hold on;

% Plot the retrieval uncertainty of the radius at cloud top
errorbar(GN_outputs.re_profile(1), GN_outputs.tau_vector(1), sqrt(GN_outputs.posterior_cov(1,1)),...
    'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the radius at cloud bottom
errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(2,2)),...
    'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the optical depth
errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(3,3)),...
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

xl0 = xline(tblut_retrieval.minRe,':',...
    ['Two-Band Look-up Table $r_{e} = $',num2str(round(tblut_retrieval.minRe, 1)), '$\mu m$'], 'Fontsize',22,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(2,:));
xl0.LabelVerticalAlignment = 'top';
xl0.LabelHorizontalAlignment = 'left';

% Plot the emit optical depth TBLUT retrieval as a constant horizontal line
yl0 = yline(tblut_retrieval.minTau,':',...
    ['Two-Band Look-up Table $\tau_{c} = $',num2str(round(tblut_retrieval.minTau, 1))], 'Fontsize',22,...
    'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(2,:));
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'right';


% compute the LWP estimate using the TBLUT retrieval
rho_liquid_water = 10^6;        % g/m^3

lwp_emit_tblut = (2*rho_liquid_water*(tblut_retrieval.minRe/1e6) * tblut_retrieval.minTau)/3; % g/m^2

% grab the hypersepctral retrieval estimate of LWP
retrieved_LWP = GN_outputs.LWP;        % g/m^2

% Print this information on the figure

dim = [.137 .35 .3 .3];
str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_emit_tblut,1)),' $g/m^{2}$', newline,...
    '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');
set(gcf,'Position',[0 0 1200 630])


% Plot the MODIS measured above cloud column water vapor





% plot the retrieved column water vapor if it was retireved
if size(GN_outputs.retrieval, 1)>3

    retrieved_CWV = GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    % Print the simulated value and the retrieved value
    str = ['$CWV_{retrieved} = \,$',num2str(round(retrieved_CWV, 2)),' $kg/m^{2}$'];

else

    % plot the assumed column water vapor used in the forward model
    % plot the HySICS simulated above cloud column water vapor
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(GN_inputs); % kg/m^2

    % print the simulated value and the foward model assumption
    str = ['$CWV_{forward \,model} = \,$',num2str(round(assumed_CWV, 2)),' $mm$'];

end


dim = [.137 .2 .3 .3];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');
set(gcf,'Position',[0 0 1200 630])


% set figure size
set(gcf,'Position',[0 0 1200 630])



end



