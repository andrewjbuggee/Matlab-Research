% Plot the upwards and downwards irradiance from the 2D monte carlo
% simulation

% By Andrew John Buggee

%%

function [] = plot_2strm_2D_monteCarlo(inputs,F_norm)

% -------------------------
% ***** Unpack Inputs *****
% -------------------------

if inputs.mie.integrate_over_size_distribution==true
    ssa = inputs.ssa_avg;
    g = inputs.g_avg;
else
    ssa = inputs.ssa;
    g = inputs.g;
end

tau_y_lower_limit = inputs.tau_y_lower_limit;
tau_y_upper_limit = inputs.tau_y_upper_limit;

% Define the albedo of the lower boundary in our model
albedo_maxTau = inputs.albedo_maxTau;


% If there are multiple layers, grab the different values
if inputs.N_layers>1
    layerRadii = inputs.layerRadii;
end




binEdges = F_norm.binEdges;


% **************************
% Check to see if there is more than 1 layer. If so, only plot the
% MonteCarlo output

if inputs.N_layers==1


    % Check to see if there is absorption

    if ssa<1

        % Next, check to see if our layer is infinitely thick, or has a finite
        % thickness

        if tau_y_upper_limit == inf

            K = sqrt((1 - ssa)*(1 - g*ssa));
            R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

            photon_fraction_up = @(tau) R_inf * exp(-K*tau);
            photon_fraction_down = @(tau) exp(-K*tau);

            % lets plot some range of tau
            tau = binEdges(1:end-1)+diff(binEdges)/2;

            C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
            C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

            figure; semilogy(photon_fraction_up(tau),tau,'Color',C1)
            hold on;
            semilogy(photon_fraction_down(tau),tau, 'Color',C2)
            semilogy(F_norm.up,tau,'Color',C1,'LineStyle',':')
            semilogy(F_norm.down,tau, 'Color',C2,'LineStyle',':')
            grid on; grid minor
            set(gca,'YDir','reverse')
            xlabel('$F/F_0$','Interpreter','latex');
            ylabel('$\tau$','Interpreter','latex');
            title({'Comparing analytical 2 stream with Monte Carlo',...
                'for an infinitely thick absorbing medium'}, 'Interpreter','latex')

            set(gcf, 'Position',[0 0 1000 630])

            % Make a text box with the parameters used
            dim = [0.65 0.65 0 0];

            texBox_str = ['$R_{\infty} = \frac{\sqrt{1 - g \tilde{\omega}}\, - \,\sqrt{1 - \tilde{\omega}}}',...
                '{\sqrt{1 - g \tilde{\omega}}\, + \,\sqrt{1 - \tilde{\omega}}}$'];
            t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
            t.Color = 'black';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'black';
            t.FitBoxToText = 'on';


            legend('$F_{\uparrow}/F_0$ 2-stream analytical','$F_{\downarrow}/F_0$ 2-stream analytical',...
                '$F_{\uparrow}/F_0$ 2D Monte Carlo','$F_{\downarrow}/F_0$ 2D Monte Carlo','Interpreter','latex',...
                'Location','best','FontSize',20)


        elseif tau_y_upper_limit>0 && tau_y_upper_limit<inf


            % K is defined in Bohren and Clothiaux (eq. 5.70)
            K = sqrt((1 - ssa)*(1 - g*ssa));
            % Define the reflectivity at the top of our layer, the photons that
            % scatter out the cloud top
            R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

            % Define the constants
            A = (R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit)/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

            B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_y_upper_limit)/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_y_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_y_upper_limit));

            photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


            photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


            % lets plot some range of tau
            tau = linspace(tau_y_lower_limit,tau_y_upper_limit,100);

            C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
            C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

            figure; plot(photon_fraction_up(tau),tau,'Color',C1)
            hold on;
            plot(photon_fraction_down(tau),tau, 'Color',C2)
            plot(F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle',':')
            plot(F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle',':')
            grid on; grid minor
            set(gca,'YDir','reverse')
            xlabel('$F/F_0$','Interpreter','latex');
            ylabel('$\tau$','Interpreter','latex');
            title({'Comparing analytical 2 stream with Monte Carlo',...
                'for an absorbing medium of finite thickness'}, 'Interpreter','latex')

            set(gcf, 'Position',[0 0 1000 630])


            dim = [0.7 0.85 0 0];

            texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
                ['N layers = ', num2str(inputs.N_layers)],...
                ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
                ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
                ['$\tilde{\omega}$ = ', num2str(round(inputs.ssa,3))], ...
                ['$g$ = ', num2str(round(inputs.g,3))],...
                ['$r$ = ', num2str(inputs.mie.radius(1)), ' $\mu m$'],...
                ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
                ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
            t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
            t.Color = 'black';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'black';
            t.FitBoxToText = 'on';

            legend('$F_{\uparrow}/F_0$ 2-stream analytical','$F_{\downarrow}/F_0$ 2-stream analytical',...
                '$F_{\uparrow}/F_0$ 2D Monte Carlo','$F_{\downarrow}/F_0$ 2D Monte Carlo','Interpreter','latex',...
                'Location','best','FontSize',20)


        end


    elseif ssa==1

        % Next, check to see if our layer is infinitely thick, or has a finite
        % thickness

        if tau_y_upper_limit == inf


            error([newline,'I dont know what to do with a layer of infinte thickness and conservative scattering.',newline])

        elseif tau_y_upper_limit>0 && tau_y_upper_limit<inf

            % For a non-absorbing layer of finite thickness, the analytical
            % solutions to the two stream radiative transfer equations are...

            R_atTop = (tau_y_upper_limit * (1 - g)/2)/(1 + tau_y_upper_limit*(1 - g)/2);

            T_atBottom = 1/(1 + tau_y_upper_limit*(1 - g)/2);

            % lets plot some range of tau
            tau = linspace(tau_y_lower_limit,tau_y_upper_limit,100);

            C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
            C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

            figure; xline(R_atTop,'Color',C1,'LineWidth',4)
            hold on;
            xline(T_atBottom,'Color',C2,'LineWidth',4)
            plot(F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle','--')
            plot(F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle','--')
            grid on; grid minor
            xlabel('$F/F_0$','Interpreter','latex');
            ylabel('$\tau$','Interpreter','latex');
            set(gca,'YDir','reverse')
            title({'Comparing analytical 2 stream with Monte Carlo',...
                'Conservative scattering for a finite layer'},...
                'Interpreter','latex')
            set(gcf, 'Position',[0 0 1000 630])

            % Make a text box with the parameters used
            dim = [0.41 0.85 0 0];

            texBox_str = {'$R(\tau = 0) = \frac{\bar{\tau}(1 - g)/2}{1 + \bar{\tau}(1 - g)/2}$',...
                '$T(\tau = \bar{\tau}) = \frac{1}{1 + \bar{\tau}(1 - g)/2}$'};
            t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
            t.Color = 'black';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'black';
            t.FitBoxToText = 'on';

            legend('$R_{\infty}$ 2-stream analytical','$T(\bar{\tau})$ 2-stream analytical',...
                '$F_{\uparrow}/F_0$ 2D Monte Carlo','$F_{\downarrow}/F_0$ 2D Monte Carlo','Interpreter','latex',...
                'Location','best','FontSize',20)


        end




    end



else


    % There is more than 1 layer in our model!! So we will ignore the
    % analytical calculation for now

    C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
    C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];
    


    % -------------------------------------------------------------------------
    % Make 2 subplots showing the results from the Monte Carlo simulation
    % and the radius value for each layer
    figure; s1 = subplot(1,2,1);
    plot(F_norm.up,binEdges(1:end-1)+diff(binEdges)/2,'Color',C1,'LineStyle',':')
    hold on;
    plot(F_norm.down,binEdges(1:end-1)+diff(binEdges)/2, 'Color',C2,'LineStyle',':')
    grid on; grid minor
    set(gca,'YDir','reverse')
    xlabel('$F/F_0$','Interpreter','latex');
    ylabel('$\tau$','Interpreter','latex');
    title({'2 Stream 1D Monte Carlo Model'}, 'Interpreter','latex')
    set(gcf, 'Position',[0 0 1200 630])


    % Make a text box with the parameters used
    dim = [0.3 0.55 0 0];

    texBox_str = {['$N_{photons}^{total} = 10^{', num2str(log10(inputs.N_photons)),'}$'],...
        ['N layers = ', num2str(inputs.N_layers)],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$\mu_0$ = ',num2str(round(cosd(inputs.solar_zenith_angle),2))],...
        ['$r_{top}$ = ',num2str(round(inputs.layerRadii(1))), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(round(inputs.layerRadii(end))), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_y_upper_limit)],...
        ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
    t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
    t.Color = 'black';
    t.FontSize = 25;
    t.FontWeight = 'bold';
    t.EdgeColor = 'black';
    t.FitBoxToText = 'on';


    % Print the legend last so the algorithm doesn't place it in front
    % of the text box
    legend('$F_{\uparrow}/F_0$ 2D Monte Carlo','$F_{\downarrow}/F_0$ 2D Monte Carlo',...
        'Interpreter','latex','Location','best','FontSize',20)


    % plot the values of each layer's particle radius at the center of each
    % layer
    d_tau = (tau_y_upper_limit - tau_y_lower_limit)/inputs.N_layers;
    tau_vector = (tau_y_lower_limit+(d_tau/2)):d_tau:(tau_y_upper_limit-(d_tau/2));

    s2 = subplot(1,2,2);
    plot(inputs.layerRadii,tau_vector,'.','MarkerSize',20)
    grid on; grid minor
    set(gca,'YDir','reverse')
    ylabel('$\tau$','Interpreter','latex');
    xlabel('$r$ $(\mu m)$', 'Interpreter','latex')
    title({'Radius of each layer'}, 'Interpreter','latex')

    % Now shift the positions of the subfigures
    s1.Position = s1.Position.*[1 1 1.5 1];

    s2.Position = s2.Position.*[1.25 1 0.5 1];


end



end