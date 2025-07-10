%% Plot the HySICS simualted droplet profile and retrieval

% INPUTS:

% (4) pixel_2Plot: a string that can either be 'first', 'median', or
% 'last'. This tells the function which pixel to use when plotting the
% MODIS results

% By Andrew John Buggee

%%

function plot_HySICS_retrieved_and_simualted_vertProf(hysics, GN_outputs, GN_inputs, tblut_retrieval)


% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
C = mySavedColors(1:2, 'fixed');


figure;



 title('Comparison between simulated and retrieved profile', 'Interpreter','latex',...
        'FontSize', 26)



% create a droplet profile from simulated measurement inputs
re_sim = create_droplet_profile2([hysics.inputs.RT.r_top, hysics.inputs.RT.r_bot],...
        GN_outputs.tau_vector, 'optical_depth', GN_inputs.model.profile.type);

tau_sim = linspace(0, hysics.inputs.RT.tau_c, 100);


% Plot the simulated measurement
plot(re_sim, tau_sim, 'Color', C(1,:),'LineStyle','-', 'LineWidth',3)

hold on

% plot two markers at cloud bottom and cloud top
plot([re_sim(1), re_sim(end)], [tau_sim(1), tau_sim(end)], '.', 'MarkerSize', 23,...
    'Color', C(1,:))

% Plot the Gauss-Newton Retrieval
plot(GN_outputs.re_profile, GN_outputs.tau_vector, 'Color', C(2,:),'LineStyle',':', 'LineWidth',3)

hold on

% plot two markers at cloud bottom and cloud top
plot([GN_outputs.re_profile(1), GN_outputs.re_profile(end)], [GN_outputs.tau_vector(1), GN_outputs.tau_vector(end)],...
    '.', 'MarkerSize', 23, 'Color', C(2,:))


% Plot the retrieval uncertainty of the radius at cloud top
errorbar(GN_outputs.re_profile(1), GN_outputs.tau_vector(1), sqrt(GN_outputs.posterior_cov(1,1)),...
    'horizontal', 'Color', C(2,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the radius at cloud bottom
errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(2,2)),...
    'horizontal', 'Color', C(2,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the optical depth
errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(3,3)),...
    'vertical', 'Color', C(2,:), 'markersize', 20, 'Linewidth', 2)



% flip the axes
set(gca, 'YDir','reverse')


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



% Plot the modis droplet estimate as a constant vertical line


xl0 = xline(tblut_retrieval.minRe,':',...
    ['TBLUT $r_{2.1} = $',num2str(round(tblut_retrieval.minRe, 2)), '$\mu m$'], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color','black');
xl0.LabelVerticalAlignment = 'bottom';


% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(tblut_retrieval.minTau,':',...
    ['TBLUT $$\tau_{2.1} = $$',num2str(round(tblut_retrieval.minTau, 2))], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color', 'black');

yl0.LabelVerticalAlignment = 'bottom';
yl0.LabelHorizontalAlignment = 'left';



% plot the HySICS simulated above cloud column water vapor
hysics_CWV_2plot = aboveCloud_CWV_simulated_hysics_spectra(hysics.inputs); % kg/m^2

% plot the retrieved column water vapor
retrieved_CWV = GN_outputs.retrieval(end, end);        % kg/m^2 (mm)



dim = [.137 .35 .3 .3];
str = ['$CWV_{simulated} = \,$',num2str(round(hysics_CWV_2plot,2)),' $kg/m^{2}$', newline,...
    '$CWV_{retrieved} = \,$',num2str(round(retrieved_CWV, 2)),' $kg/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');
set(gcf,'Position',[0 0 1200 630])

% Create a Legend with only the two black curves
%legend('Vocals Rex In-situ Measurement', 'Desired Retrieval Profile', 'Interpreter','latex', 'Location','best')
legend({'Simulated $r_e(\tau)$', '', 'Retrieved $r_e(\tau)$'}, 'Interpreter','latex', 'Location','northwest', 'FontSize', 25)



% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, abs(GN_inputs.RT.z(end) - GN_inputs.RT.z(1))])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left

grid on; grid minor

ylim([-0.1, 1.2*GN_inputs.RT.tau_c])

end