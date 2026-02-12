%% Plot the EMIT retrieved droplet profile and optical depth for a single pixel
% Also plot the coincident retrievals from MODIS, AIRS, and AMSR-E

% INPUTS:



% By Andrew John Buggee

%%

function figure_handle = plot_EMIT_retrieved_vertProf_with_MODIS_AIRS_AMSR_perPixel(GN_outputs, GN_inputs,...
    modis, airs, amsr, pixel_num)



if isempty(airs) == false

    % define where the US standard atm's are located

    % Determine which computer you're using
    which_computer = whatComputer();

    % Load simulated measurements
    if strcmp(which_computer,'anbu8374')==true

        % ------ Folders on my Mac Desktop --------
        % -----------------------------------------

        atm_data_directory = ['/Users/andrewbuggee/Documents/libRadtran-2.0.6/data/atmmod/'];

    elseif strcmp(which_computer,'andrewbuggee')==true

        % ------ Folders on my Macbook --------
        % -------------------------------------
        atm_data_directory = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/'];

    end


end








C = mySavedColors(1:2, 'fixed');

% Plot the Gauss-Newton Retrieval

figure_handle = figure;

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
if isfield(GN_outputs, 'posterior_cov')==true


    errorbar(GN_outputs.re_profile(1), GN_outputs.tau_vector(1), sqrt(GN_outputs.posterior_cov(1,1)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the radius at cloud bottom
    errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(2,2)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the optical depth
    errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov(3,3)),...
        'vertical', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)


elseif isfield(GN_outputs, 'posterior_cov_log')==true


    errorbar(GN_outputs.re_profile(1), GN_outputs.tau_vector(1), sqrt(GN_outputs.posterior_cov_log(1,1)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the radius at cloud bottom
    errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov_log(2,2)),...
        'horizontal', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)

    % Plot the retrieval uncertainty of the optical depth
    errorbar(GN_outputs.re_profile(end), GN_outputs.tau_vector(end), sqrt(GN_outputs.posterior_cov_log(3,3)),...
        'vertical', 'Color',mySavedColors(1,'fixed'), 'markersize', 20, 'Linewidth', 2)


end



% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.00857142857142857 0.884832451499119 0.051 0.0777777777777779],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.00857142857142857 0.0262081128747797 0.051 0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');




% ----------------------------------------------------------------
% ------------------ plot MODIS TBLUT retrieval ------------------
% ----------------------------------------------------------------


% Plot the emit TBLUT droplet estimate as a constant vertical line

xl0 = xline(modis.cloud.effRadius17(pixel_num),':',...
    ['MODIS retreived $r_{e} = $',num2str(round(modis.cloud.effRadius17(pixel_num), 1)), '$\mu m$'], 'Fontsize',22,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(2,:));
xl0.LabelVerticalAlignment = 'top';
xl0.LabelHorizontalAlignment = 'left';

% Plot the emit optical depth TBLUT retrieval as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(pixel_num),':',...
    ['MODIS retrieved $\tau_{c} = $',num2str(round(modis.cloud.optThickness17(pixel_num), 1))], 'Fontsize',22,...
    'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(2,:));
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'right';




% Plot the emit TBLUT droplet estimate as a constant vertical line

% xl0 = xline(tblut_retrieval.minRe,':',...
%     ['Two-Band Look-up Table $r_{e} = $',num2str(round(tblut_retrieval.minRe, 1)), '$\mu m$'], 'Fontsize',22,...
%     'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(2,:));
% xl0.LabelVerticalAlignment = 'top';
% xl0.LabelHorizontalAlignment = 'left';
%
% % Plot the emit optical depth TBLUT retrieval as a constant horizontal line
% yl0 = yline(tblut_retrieval.minTau,':',...
%     ['Two-Band Look-up Table $\tau_{c} = $',num2str(round(tblut_retrieval.minTau, 1))], 'Fontsize',22,...
%     'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(2,:));
% yl0.LabelVerticalAlignment = 'top';
% yl0.LabelHorizontalAlignment = 'right';


% % compute the LWP estimate using the TBLUT retrieval
% rho_liquid_water = 10^6;        % g/m^3
%
% lwp_emit_tblut = (2*rho_liquid_water*(tblut_retrieval.minRe/1e6) * tblut_retrieval.minTau)/3; % g/m^2





% grab the hypersepctral retrieval estimate of LWP
retrieved_LWP = GN_outputs.LWP;        % g/m^2

% ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3

lwp_modis_WH = 5/9 * rho_h2o * modis.cloud.optThickness17(pixel_num) *...
    (modis.cloud.effRadius17(pixel_num) * 1e-6);

% Print this information on the figure
dim = [0.137 0.727962757628641 0.481571428571429 0.159691563359013];

if isempty(amsr) == false & isnan(amsr.cloud.LiquidWaterPath(pixel_num))

    disp([newline, 'AMSR-E data isnt valid at this paxiel: NaN', newline])

    str = ['$LWP_{MODIS} = \,$',num2str(round(modis.cloud.lwp(pixel_num),1)),' $g/m^{2}$', newline,...
        '$LWP_{MODIS-WH} = \,$',num2str(round(lwp_modis_WH,1)),' $g/m^{2}$', newline,...
        '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];


elseif isempty(amsr) == false

    str = ['$LWP_{MODIS} = \,$',num2str(round(modis.cloud.lwp(pixel_num),1)),' $g/m^{2}$', newline,...
        '$LWP_{MODIS-WH} = \,$',num2str(round(lwp_modis_WH,1)),' $g/m^{2}$', newline,...
        '$LWP_{AMSR} = \,$',num2str(round(amsr.cloud.LiquidWaterPath(pixel_num),1)),' $g/m^{2}$', newline,...
        '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];

else

    str = ['$LWP_{MODIS} = \,$',num2str(round(modis.cloud.lwp(pixel_num),1)),' $g/m^{2}$', newline,...
        '$LWP_{MODIS-WH} = \,$',num2str(round(lwp_modis_WH,1)),' $g/m^{2}$', newline,...
        '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$'];


end


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');


% Plot the MODIS measured above cloud column water vapor





% plot the retrieved column water vapor if it was retireved
if size(GN_outputs.retrieval, 1)>3

    retrieved_CWV = GN_outputs.retrieval(end, end);        % kg/m^2 (mm)



    % Print the simulated value and the retrieved value
    % MODIS reports in cm
    % My retrieval is is kg/m^2 which is equivelant to mm

    if isempty(airs) == true || (isempty(airs) == false & isnan(airs.H2O.totCol_Std(pixel_num)))

        disp([newline, 'AMSR-E data isnt valid at this paxiel: NaN', newline])

        str = ['$ACPW_{MODIS} = \,$',num2str(round(modis.vapor.col_nir(pixel_num) * 10, 1)),' $mm$', newline,...
            '$ACPW_{Hyperspectral} = \,$',num2str(round(retrieved_CWV, 1)),' $mm$'];

    elseif isempty(airs) == false


        str = ['$ACPW_{MODIS} = \,$',num2str(round(modis.vapor.col_nir(pixel_num) * 10, 1)),' $mm$', newline,...
            '$ACPW_{AIRS} = \,$',num2str(round(airs.H2O.acpw_using_assumed_CTH(pixel_num), 1)),' $mm$', newline,...
            '$ACPW_{Hyperspectral} = \,$',num2str(round(retrieved_CWV, 1)),' $mm$'];

    end

else

    % plot the assumed column water vapor used in the forward model
    % plot the HySICS simulated above cloud column water vapor
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(GN_inputs); % kg/m^2

    % print the simulated value and the foward model assumption
    str = ['$ACPW_{forward \,model} = \,$',num2str(round(assumed_CWV, 2)),' $mm$'];

end


dim = [0.137 0.559773612316743 0.475857142857143 0.118621449411651];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');


% set figure size
set(gcf,'Position',[0 0 700 810])



end



