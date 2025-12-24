% Plot single instance of in-situ measurement and retreival

% By Andrew John Buggee

function plot_retrieved_prof_with_inSitu_paper2(mat_file_path, mat_file_name)


% Load the data from the file
ds = load([mat_file_path, mat_file_name]);




ln_wdth = 1;
mkr_sz = 20;

C = mySavedColors(61:68, 'fixed');
C_idx = 5;


% Plot the in-situ profile

figure;

title('In-situ vs. Retrieved droplet profile', 'Interpreter','latex',...
    'FontSize', 26)

if ds.GN_inputs.measurement.tau(end) ~= ds.GN_inputs.measurement.tau_c

    % check the length
    if length(ds.GN_inputs.measurement.tau) ~= length(ds.GN_inputs.measurement.re_prof)

        warning([newline, 'Tau vector length doesnt match radius profile length', newline])

        if ds.GN_inputs.measurement.tau(length(ds.GN_inputs.measurement.re_prof)) == ds.GN_inputs.measurement.tau_c

            plot(ds.GN_inputs.measurement.re_prof, ds.GN_inputs.measurement.tau(1:length(ds.GN_inputs.measurement.re_prof)),...
                'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
                'Color', 'k', 'MarkerFaceColor', 'k')


        end

    end

else


    plot(ds.GN_inputs.measurement.re_prof, ds.GN_inputs.measurement.tau,...
        'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
        'Color', 'k', 'MarkerFaceColor', 'k')

end

% flip y-axis and provide axes labels
set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex', 'FontSize', 35)
grid on; grid minor; hold on;


% plot the retrieved profile
plot(ds.GN_outputs.re_profile, ds.GN_outputs.tau_vector,...
    'LineStyle','-', 'LineWidth',ln_wdth+2,'Color', C(C_idx,:))

% Plot the retrieval uncertainty of the radius at cloud top
errorbar(ds.GN_outputs.re_profile(1), ds.GN_outputs.tau_vector(1), sqrt(ds.GN_outputs.posterior_cov_log(1,1)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the radius at cloud bottom
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(2,2)),...
    'horizontal', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)

% Plot the retrieval uncertainty of the optical depth
errorbar(ds.GN_outputs.re_profile(end), ds.GN_outputs.tau_vector(end), sqrt(ds.GN_outputs.posterior_cov_log(3,3)),...
    'vertical', 'Color', C(C_idx,:), 'markersize', 20, 'Linewidth', 2)



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



% Plot the TBLUT droplet estimate as a constant vertical line

xl0 = xline(ds.tblut_retrieval.minRe,':',...
    ['TBLUT $r_{e} = $',num2str(round(ds.tblut_retrieval.minRe, 1)), '$\mu m$'], 'Fontsize',24,...
    'FontWeight', 'bold', 'Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
xl0.LabelVerticalAlignment = 'bottom';
xl0.LabelHorizontalAlignment = 'left';

% Plot the optical depth TBLUT retrieval as a constant horizontal line
yl0 = yline(ds.tblut_retrieval.minTau,':',...
    ['TBLUT $\tau_{c} = $',num2str(round(ds.tblut_retrieval.minTau, 1))], 'Fontsize',24,...
    'FontWeight', 'bold','Interpreter','latex','LineWidth',3,'Color', C(C_idx,:));
yl0.LabelVerticalAlignment = 'top';
yl0.LabelHorizontalAlignment = 'right';


% compute the LWP estimate using the TBLUT retrieval
rho_liquid_water = 10^6;        % g/m^3

lwp_tblut = (2*rho_liquid_water*(ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

% grab the hypersepctral retrieval estimate of LWP
retrieved_LWP = ds.GN_outputs.LWP;        % g/m^2

% What is the true LWP
LWP_true = ds.GN_inputs.measurement.lwp;   % g/m^2

% Print this information on the figure

dim = [0.141166666666667 0.71690449790349 0.293506310780843 0.157698676699684];
str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_tblut,1)),' $g/m^{2}$', newline,...
    '$LWP_{hyperspectral} = \,$',num2str(round(retrieved_LWP,1)),' $g/m^{2}$', newline...
    '$LWP_{true} = \,$',num2str(round(LWP_true,1)),' $g/m^{2}$'];

annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');




% plot the retrieved column water vapor if it was retireved
if size(ds.GN_outputs.retrieval, 1)>3

    acpw_true = ds.GN_inputs.measurement.actpw;             % kg/m^2 (mm)
    retrieved_CWV = ds.GN_outputs.retrieval(end, end);        % kg/m^2 (mm)

    % Print the simulated value and the retrieved value
    str = ['$acpw_{true} = \,$',num2str(round(acpw_true, 2)),' $mm$', newline,...
        '$acpw_{retrieved} = \,$',num2str(round(retrieved_CWV, 2)),' $mm$'];

else

    % plot the assumed column water vapor used in the forward model
    % plot the HySICS simulated above cloud column water vapor
    assumed_CWV = aboveCloud_CWV_simulated_hysics_spectra(ds.GN_inputs); % kg/m^2

    % print the simulated value and the foward model assumption
    str = ['$acpw_{forward \,model} = \,$',num2str(round(assumed_CWV, 2)),' $mm$', newline,...
        '$acpw_{MODIS} = \,$',num2str(modis_retrieved_aboveCloud_CWV),' $mm$'];

end


dim = [0.141166666666667 0.571252065991597 0.238002173105876 0.103351108611576];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');

title(['Retrievel using all 636 HySICS spectral channels'],...
    'Fontsize', 25, 'Interpreter', 'latex');

legend('In-Situ', 'Hyperspectral Retrieval', 'Interpreter','latex', 'Location','best', 'FontSize', 20,...
            'Color', 'white', 'TextColor', 'k')


% set figure size
set(gcf,'Position',[0 0 1200 630])


end