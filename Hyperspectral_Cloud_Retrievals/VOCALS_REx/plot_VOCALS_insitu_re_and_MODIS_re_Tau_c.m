% Plot vocals-rex in-situ meausrement of effective droplet radius versus
% optical depth. Overlap the MODIS retrieval for effective radius and for
% cloud optical depth

function [] = plot_VOCALS_insitu_re_and_MODIS_re_Tau_c(modis,vocalsRex)


% ------------------------------------------------------------------
% ------------------ Compute optical depth -------------------------
% ------------------------------------------------------------------
% optical depth is defined to be 0 at cloud top and increasing towards
% cloud bottom
vector_length = length(vocalsRex.altitude);
tau = zeros(1,vector_length-1);
for ii = 1:vector_length-1
    
    % we have to convert Nc and re to have the same units as the alitude,
    % which is in meters
    re_meters = vocalsRex.re(vector_length-ii:vector_length)./1e6;                                   % meters
    total_Nc_meters = vocalsRex.total_Nc(vector_length-ii:vector_length).*1e6;                           % #/m^3
    altitude = vocalsRex.altitude(end) -  vocalsRex.altitude(vector_length-ii:vector_length);
    % we need to flip these vectors so we start integratin at the cloud
    % top!
    tau(ii) = 2*pi* trapz(fliplr(altitude), fliplr(re_meters.^2 .* total_Nc_meters));

end

% add a zero at the begining!
tau = [0,tau];


% ------------------------------------------------------------------
% ---------------- Compute Liquid Water Path -----------------------
% ------------------------------------------------------------------

LWP = trapz(vocalsRex.altitude, vocalsRex.lwc);                 % grams of water/m^2



% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];
nice_orange = [0.8500, 0.3250, 0.0980];

figure;
plot(fliplr(vocalsRex.re),tau,'-o','Color','black', 'MarkerSize',10,...
    'MarkerFaceColor','black','LineWidth',1); 
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35); 
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
title('Comparison between in-situ and MODIS retrieved $r_e$', 'Interpreter','latex')
grid on; grid minor; hold on; 

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

% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.029,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.029,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the modis droplet estimate as a constant vertical line
xl0 = xline(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist),':',['MODIS $$r_{2.1} = $$',num2str(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist)), '$$\mu m$$'], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist),':',['MODIS $$\tau_{2.1} = $$',num2str(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist))], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';


% Let's compute the mean number concentration within this cloud and print
% it on our plot

mean_Nc = mean(vocalsRex.total_Nc);

dim = [.2 .5 .3 .3];
str = ['$$< N_c > = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$',newline,'$$LWP = $$',num2str(round(LWP,1)),' $$g/m^{2}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18,'FontWeight','bold');
set(gcf,'Position',[0 0 1100 630])

% Create a Legend with only the two black curves
%legend('Vocals Rex In-situ Measurement', 'Desired Retrieval Profile', 'Interpreter','latex', 'Location','best')
legend('Vocals Rex In-situ Measurement', 'Interpreter','latex', 'Location','best')


end
