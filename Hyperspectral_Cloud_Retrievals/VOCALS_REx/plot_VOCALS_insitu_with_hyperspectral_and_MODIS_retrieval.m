function [] = plot_VOCALS_insitu_with_hyperspectral_and_MODIS_retrieval(modis,vocalsRex, GN_outputs)


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

% Plot the vocals rex data as markers, but add a line through each to guide
% the eye
plot(fliplr(vocalsRex.re),tau,'Marker','.','MarkerSize',35, 'LineStyle','-', 'LineWidth',1);
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex', 'FontSize', 35); 
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
title('Comparison of Gauss-Newton and MODIS retrieved $r_e$ with in-situ', 'Interpreter','latex')
grid on; grid minor; hold on; 

% Plot the Gauss-Newton Retreval

hyperspectral_retrieval_linewidth = 6;
hyperspectral_retrieval_color = mySavedColors(1,'fixed');                        % Bright pink

plot(GN_outputs.re_profile, GN_outputs.tau_vector, 'Color',hyperspectral_retrieval_color,'LineStyle',...
    ':', 'LineWidth',hyperspectral_retrieval_linewidth)

% Plot the modis droplet estimate as a constant vertical line
xl0 = xline(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist),':',...
    ['MODIS $$r_{2.1} = $$',num2str(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist)), '$$\mu m$$'],...
    'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'middle';
xl0.LabelHorizontalAlignment = 'right';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist),':',...
    ['MODIS $$\tau_{2.1} = $$',num2str(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist))],...
    'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'left';


% Plot the our GAuss-Newton retireval for optical depth as a constant horizontal line
yl1 = yline(GN_outputs.retrieval(3,end),':',...
    ['Hyperspectral retrieved $$\tau_{c} = $$',num2str(round(GN_outputs.retrieval(3,end),2))],...
    'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',hyperspectral_retrieval_color,'LabelVerticalAlignment','bottom',...
    'LabelHorizontalAlignment','center');

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




% Let's compute the mean number concentration within this cloud and print
% it on our plot
mean_Nc = mean(vocalsRex.total_Nc);

dim = [.2 .5 .3 .3];
str = ['$$< N_c > = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$',newline,'$$LWP = $$',num2str(round(LWP,1)),' $$g/m^{2}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18,'FontWeight','bold');

% set figure size
set(gcf,'Position',[0 0 1300 730])

% Create a Legend with only the two black curves
legend('Vocals Rex In-situ Measurement', 'Retrieved Profile with first 7 channels of MODIS', 'Interpreter','latex', 'Location','best')



end