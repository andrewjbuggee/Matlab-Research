%% Plot the VOCALS-Rex in-situ data along with a MODIS retrieved droplet size and optical depth for a single pixel

% INPUTS:

% (4) pixel_2Plot: a string that can either be 'first', 'median', or
% 'last'. This tells the function which pixel to use when plotting the
% MODIS results

% By Andrew John Buggee

%%

function plot_vocalsRex_with_MODIS_retrieved_re_and_vertProf_retrieval(vocalsRex, modis, modisInputs, GN_outputs, GN_inputs, pixel_2Plot)


% ------------------------------------------------------------------
% ---------------- Compute Liquid Water Path -----------------------
% ------------------------------------------------------------------

LWP_vocals = trapz(vocalsRex.altitude, vocalsRex.lwc);                 % grams of water/m^2

% ----- Compute the CDP uncertainty -----
re_uncertainty = cloud_droplet_probe_uncertainty_estimate(vocalsRex.re);



% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];
nice_orange = [0.8500, 0.3250, 0.0980];

figure;
errorbar(flipud(vocalsRex.re), vocalsRex.tau', flipud(re_uncertainty), 'horizontal','-o','Color','black', 'MarkerSize',10,...
    'MarkerFaceColor','black','LineWidth',1);
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
grid on; grid minor; hold on;

if modisInputs.flags.useAdvection==true

    title('With Advection - Comparison between in-situ, MODIS retrieval and vertical retrieval', 'Interpreter','latex',...
        'FontSize', 26)

else

    title('Without Advection - Comparison between in-situ, MODIS retrieval and vertical retrieval', 'Interpreter','latex',...
        'FontSize', 26)

end

legend_str = cell(1, 3*size(GN_outputs.re_profile, 2));
idx_step = 0;

rms_values  = zeros(1, size(GN_outputs.re_profile,2));

% Plot the Gauss-Newton Retrieval
for pp = 1:size(GN_outputs.re_profile,2)

    plot(GN_outputs.re_profile(:,pp), GN_outputs.tau_vector(:,pp), 'Color',mySavedColors(pp,'fixed'),'LineStyle',':', 'LineWidth',3)


    % Plot the retrieval uncertainty of the radius at cloud top
    errorbar(GN_outputs.re_profile(1,pp), GN_outputs.tau_vector(1,pp), sqrt(GN_outputs.posterior_cov(1,1,pp)),...
        'horizontal', 'Color',mySavedColors(pp,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the radius at cloud bottom
    errorbar(GN_outputs.re_profile(end,pp), GN_outputs.tau_vector(end,pp), sqrt(GN_outputs.posterior_cov(2,2,pp)),...
        'horizontal', 'Color',mySavedColors(pp,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the optical depth
    errorbar(GN_outputs.re_profile(end,pp), GN_outputs.tau_vector(end,pp), sqrt(GN_outputs.posterior_cov(3,3,pp)),...
        'vertical', 'Color',mySavedColors(pp,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Store the rms residual
    rms_values(pp) = GN_outputs.rms_residual{pp}(end);

    % create legend string for each index
%     legend_str{3*pp - 2} = ['Retrieved Profile - idx = ', num2str(pp), ' - rms residual = ', ...
%         num2str(rms_values(pp))];
    legend_str{4*pp - 3} = ['idx = ', num2str(pp), ' - rms residual = ', ...
            num2str(round(rms_values(pp),5))];

    % create empty legends for the error bar entries
    idx_step = idx_step+1;
    legend_str{pp+idx_step} = '';

    idx_step = idx_step+1;
    legend_str{pp+idx_step} = '';

    idx_step = idx_step+1;
    legend_str{pp+idx_step} = '';


end


% Fit a curve to the in-situ data to show the capability we are interested
% in devloping

% curve_fit_linewidth = 6;
% curve_fit_color = mySavedColors(1,'fixed');                        % Bright pink
%
% f = fit(tau', fliplr(double(vocalsRex.re))', 'smoothingspline','SmoothingParam',0.9);
% hold on;
% plot(f(tau),tau,'Color',curve_fit_color,'LineStyle',':', 'LineWidth',curve_fit_linewidth);

% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, vocalsRex.altitude(end) - vocalsRex.altitude(1)])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left
ylim([-0.25, vocalsRex.tau(end)*1.15])

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


xl0 = xline(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(pixel_2Plot)),':',...
    ['MODIS $$r_{2.1} = $$',num2str(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(pixel_2Plot))), '$$\mu m$$'], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(pixel_2Plot)),':',...
    ['MODIS $$\tau_{2.1} = $$',num2str(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(pixel_2Plot)))], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_orange);



yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';



% grab the MODIS LWP to plot
modis_lwp_2plot = modis.cloud.lwp(vocalsRex.modisIndex_minDist(pixel_2Plot)); % g/m^2

% grab the retrieved LWP to plot
retrieved_LWP = GN_outputs.LWP(pixel_2Plot);        % g/m^2

% Let's compute the mean number concentration within this cloud and print
% it on our plot

mean_Nc = mean(vocalsRex.Nc);

dim = [.137 .35 .3 .3];
str = ['$$< N_c >_{in-situ} = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$',newline,...
    '$$LWP_{in-situ} = \,$$',num2str(round(LWP_vocals,1)),' $$g/m^{2}$$',newline,...
    '$LWP_{MODIS} = \,$',num2str(round(modis_lwp_2plot,1)),' $g/m^{2}$', newline,...
    '$LWP_{retrieved} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',20,'FontWeight','bold');
set(gcf,'Position',[0 0 1200 630])

% Create a Legend with only the two black curves
%legend('Vocals Rex In-situ Measurement', 'Desired Retrieval Profile', 'Interpreter','latex', 'Location','best')
legend(['Vocals Rex In-situ Measurement', legend_str], 'Interpreter','latex', 'Location','northwest', 'FontSize', 20)


% Include a text box stating the percentage of the TBLUT guess that is used
% as the standard deviation
% Create textbox
annotation('textbox',[0.6893 0.0142 0.27 0.04126],...
    'String',{['$r_{top}$ = ', num2str(sqrt(GN_inputs.model.covariance(1,1,1))/GN_inputs.model.apriori(1,1)*100), '$\%$',...
    ' $ \; \; \; r_{bot}$ = ', num2str(sqrt(GN_inputs.model.covariance(2,2,1))/GN_inputs.model.apriori(1,2)*100), '$\%$',...
    ' $\; \; \; \tau_c$ = ', num2str(sqrt(GN_inputs.model.covariance(3,3,1))/GN_inputs.model.apriori(1,3)*100), '$\%$']},...
    'FitBoxToText','on',...
    'Interpreter','latex',...
    'LineWidth', 2,...
    'FontSize',20,'FontWeight','bold');


end