% Plot single instance of in-situ measurement and retreival

% By Andrew John Buggee

function fig1 = plot_retrieved_prof_with_inSitu_paper2(mat_file_path, mat_file_name)





% Determine which computer you're using
which_computer = whatComputer();

% Load simulated measurements
if strcmp(which_computer,'anbu8374')==true

    % -----------------------------------------
    % ------ Folders on my Mac Desktop --------
    % -----------------------------------------

elseif strcmp(which_computer,'andrewbuggee')==true

    % -------------------------------------
    % ------ Folders on my Macbook --------
    % -------------------------------------

    % mie folder location
    mie_folder_path = '/Users/andrewbuggee/Documents/libRadtran-2.0.6/Mie_Calculations/';


elseif strcmp(which_computer,'curc')==true


    % ------------------------------------------------
    % ------ Folders on the CU Super Computer --------
    % ------------------------------------------------


end






% Load the data from the file
ds = load([mat_file_path, mat_file_name]);




ln_wdth = 1;
mkr_sz = 20;

C = mySavedColors(61:68, 'fixed');
C_idx = 5;


% Plot the in-situ profile

fig1 = figure;

% Create axes
axes1 = axes('Parent',fig1,'Position',[0.175438596491228 0.11 0.804561403508772 0.815]);
hold(axes1,'on');

title('In-situ vs. Retrieved droplet profile', 'Interpreter','latex',...
    'FontSize', 26)

if isfield(ds.GN_inputs.measurement, 'tau')

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






elseif isfield(ds.GN_inputs.measurement, 'tau_prof')==true


    % check the length
    if length(ds.GN_inputs.measurement.tau_prof) == length(ds.GN_inputs.measurement.re_prof)


        % ** tau_prof starts at the bottom and re_prof starts at cloud top
        % **

        plot(ds.GN_inputs.measurement.re_prof, sort(ds.GN_inputs.measurement.tau_prof),...
            'Marker','.','LineStyle','-', 'LineWidth',ln_wdth, 'MarkerSize', mkr_sz,...
            'Color', 'k', 'MarkerFaceColor', 'k')

    else

        error([newline, 'Tau vector length doesnt match radius profile length', newline])

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
annotation('textbox',[0.0109176753017712 0.862655122655124 0.0509999999999998 0.0777777777777779],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.0109176753017712 0.0834920634920637 0.051 0.0777777777777779],...
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
yl0.LabelVerticalAlignment = 'bottom';
yl0.LabelHorizontalAlignment = 'left';


% compute the LWP estimate using the TBLUT retrieval
con = physical_constants;
rho_h2o = con.density_h2o_liquid * 1e3;   % g/m^3

lwp_tblut = (2 * rho_h2o * (ds.tblut_retrieval.minRe/1e6) * ds.tblut_retrieval.minTau)/3; % g/m^2

% ** Compute the Wood-Hartmann LWP estimate asssuming Adiabatic **
lwp_tblut_WH = 5/9 * rho_h2o * ds.tblut_retrieval.minTau * (ds.tblut_retrieval.minRe/1e6);


% grab the hypersepctral retrieval estimate of LWP
% retrieved_LWP = ds.GN_outputs.LWP;        % g/m^2
% -------------------------------------------------------
% -------------------------------------------------------
% ** Compute new updated LWP calc ***
re_profile = create_droplet_profile2([ds.GN_outputs.retrieval(1,end), ds.GN_outputs.retrieval(2,end)],...
    ds.GN_inputs.RT.z, 'altitude', ds.GN_inputs.model.profile.type);


% define the z vector
z = linspace(ds.GN_inputs.RT.z_topBottom(2), ds.GN_inputs.RT.z_topBottom(1), length(re_profile)+1)';                 % km - altitude vector

% define the z midpoint at each layer and normalize it!
z_norm = z - z(1);
z_norm_mid = (diff(z_norm)/2 + z_norm(1:end-1));


% The radius input is defined as [r_start, r_end, r_step].
% where r_step is the interval between radii values (used only for
% vectors of radii). A 0 tells the code there is no step. Finally, the
% radius values have to be in increasing order.
ext_bulk_coeff_per_LWC = zeros(length(re_profile), 1);

for rr = 1:length(re_profile)

    mie_radius = [re_profile(rr), re_profile(rr), 0];    % microns

    size_distribution = {'gamma', ds.GN_inputs.RT.distribution_var(rr)};           % droplet distribution

    % Create a mie file
    [input_filename, output_filename] = write_mie_file('MIEV0', 'water',...
        mie_radius, 500, size_distribution, 'verbose', rr, round(re_profile(rr), 4), mie_folder_path);

    % run the mie file
    [~] = runMIE(mie_folder_path, input_filename,output_filename, which_computer);

    % Read the output of the mie file
    [mie,~,~] = readMIE(mie_folder_path, output_filename);

    ext_bulk_coeff_per_LWC(rr) = mie.Qext;       % km^-1 / (cm^3 / m^3)

end


% ** Assuming liquid water content increases linearly with depth **
z_kilometers_upper_boundary = z(2:end) - z(1);                     % kilometers - geometric depth at upper boundary of each cloud layer
dz_km = z(2) - z(1);           % kilometers

%slope = tau_c /(dz_km * sum(ext_bluk_coeff_per_LWC .* z_kilometers_midpoint ));     % g/m^3/m - slope of the lwc profile
slope = ds.GN_outputs.retrieval(3,end) /(dz_km * sum(ext_bulk_coeff_per_LWC .* z_kilometers_upper_boundary ));     % g/m^3/km - slope of the lwc profile

% solve for the linear liquid water content profile
%lwc = slope * z_kilometers_midpoint;                     % g/m^3 - grams of water per meter cubed of air
lwc = slope * z_kilometers_upper_boundary;                     % g/m^3 - grams of water per meter cubed of air


retrieved_LWP = trapz( 1e3 .* z_norm_mid, lwc);    % g/m^2


% -------------------------------------------------------
% -------------------------------------------------------



% What is the true LWP
LWP_true = ds.GN_inputs.measurement.lwp;   % g/m^2

% Print this information on the figure

dim = [0.185672925359295 0.770354933354154 0.477184217497848 0.122654526321976];
str = ['$LWP_{TBLUT} = \,$',num2str(round(lwp_tblut,1)),' $g/m^{2}$', newline,...
    '$LWP_{TBLUT-WH} = \,$',num2str(round(lwp_tblut_WH,1)),' $g/m^{2}$', newline,...
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


dim = [0.185672925359295 0.633859831990577 0.385755646069276 0.0803841955867814];


annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',25,'FontWeight','bold');

title(['Retrievel using all 636 HySICS spectral channels'],...
    'Fontsize', 25, 'Interpreter', 'latex');

legend('In-Situ', 'Hyperspectral Retrieval', 'Interpreter','latex', 'Position',...
    [0.191290035273543 0.536882716049383 0.335733019871267 0.0648148148148148], 'FontSize', 20,...
    'Color', 'white', 'TextColor', 'k')


% set figure size
set(gcf,'Position',[0 0 700 810])

box(axes1,'on');
grid(axes1,'on');
axis(axes1,'ij');
hold(axes1,'off');

% set x axis limits
xlim([min(ds.GN_inputs.measurement.re_prof) - 0.5, max(ds.GN_inputs.measurement.re_prof) + 0.5])


end