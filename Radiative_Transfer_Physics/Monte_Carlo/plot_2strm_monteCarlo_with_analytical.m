% Plot and compare the two stream monte carlo solutions to the analytical
% solution defined in Bohren and Clothiaux

% By Andrew John Buggee

%%

function [] = plot_2strm_monteCarlo_with_analytical(inputs,F_norm)

% -------------------------
% ***** Unpack Inputs *****
% -------------------------

ssa = inputs.ssa;
g = inputs.g;

tau_lower_limit = inputs.tau_lower_limit;
tau_upper_limit = inputs.tau_upper_limit;

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

        if tau_upper_limit == inf

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
            t.Color = 'white';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'white';
            t.FitBoxToText = 'on';


            legend('$F_{\uparrow}/F_0$ analytical','$F_{\downarrow}/F_0$ analytical',...
                '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
                'Location','best')


        elseif tau_upper_limit>0 && tau_upper_limit<inf


            % K is defined in Bohren and Clothiaux (eq. 5.70)
            K = sqrt((1 - ssa)*(1 - g*ssa));
            % Define the reflectivity at the top of a semi-infinite cloud
            % R_inf = F_up(tau=0)/F_0 
            % This is the fraction of light that scatters out the cloud top
            R_inf = (sqrt(1-ssa*g) - sqrt(1 - ssa))/(sqrt(1-ssa*g) + sqrt(1 - ssa));

            % Define the constants
            A = (R_inf - albedo_maxTau)*exp(-K*tau_upper_limit)/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_upper_limit));

            B = -(1 - R_inf*albedo_maxTau)*exp(K*tau_upper_limit)/...
                (R_inf*(R_inf - albedo_maxTau)*exp(-K*tau_upper_limit) - (1 - albedo_maxTau*R_inf)*exp(K*tau_upper_limit));

            photon_fraction_up = @(tau) A*exp(K*tau) + B*R_inf*exp(-K*tau);


            photon_fraction_down = @(tau) A*R_inf*exp(K*tau) + B*exp(-K*tau);


            % lets plot some range of tau
            tau = linspace(tau_lower_limit,tau_upper_limit,100);

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

            texBox_str = {[num2str(inputs.N_photons),' photons'],...
                ['N layers = ', num2str(inputs.N_layers)],...
                ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
                ['$\tilde{\omega}$ = ', num2str(inputs.ssa)], ...
                ['$g$ = ', num2str(inputs.g)],...
                ['$r$ = ', num2str(inputs.mie.radius(1)), ' $\mu m$'],...
                ['$\tau_0$ = ', num2str(inputs.tau_upper_limit)],...
                ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
            t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
            t.Color = 'white';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'white';
            t.FitBoxToText = 'on';

            legend('$F_{\uparrow}/F_0$ analytical','$F_{\downarrow}/F_0$ analytical',...
                '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
                'Location','best')


        end


    elseif ssa==1

        % Next, check to see if our layer is infinitely thick, or has a finite
        % thickness

        if tau_upper_limit == inf


            error([newline,'I dont know what to do with a layer of infinte thickness and conservative scattering.',newline])

        elseif tau_upper_limit>0 && tau_upper_limit<inf

            % For a non-absorbing layer of finite thickness, the analytical
            % solutions to the two stream radiative transfer equations are...

            R_atTop = (tau_upper_limit * (1 - g)/2)/(1 + tau_upper_limit*(1 - g)/2);

            T_atBottom = 1/(1 + tau_upper_limit*(1 - g)/2);

            % lets plot some range of tau
            tau = linspace(tau_lower_limit,tau_upper_limit,100);

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
            t.Color = 'white';
            t.FontSize = 25;
            t.FontWeight = 'bold';
            t.EdgeColor = 'white';
            t.FitBoxToText = 'on';

            legend('$R_{\infty}$ analytical','$T(\bar{\tau})$ analytical',...
                '$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo','Interpreter','latex',...
                'Location','best')


        end




    end



else


    % There is more than 1 layer in our model!! So we will ignore the
    % analytical calculation for now

    C1 = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01];
    C2 = [4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01];

    % Make 2 subplots showing the results from the Monte Carlo simulation
    % and the values used for the cloud layers
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

    texBox_str = {[num2str(inputs.N_photons),' photons'],...
        ['N layers = ', num2str(inputs.N_layers)],...
        ['$\lambda$ = ',num2str(inputs.mie.wavelength(1)), ' $nm$'],...
        ['$r_{top}$ = ',num2str(inputs.layerRadii(1)), ' $\mu m$'],...
        ['$r_{bot}$ = ',num2str(inputs.layerRadii(end)), ' $\mu m$'],...
        ['$\tau_0$ = ', num2str(inputs.tau_upper_limit)],...
        ['$A_0$ = ', num2str(inputs.albedo_maxTau)]};
    t = annotation('textbox',dim,'string',texBox_str,'Interpreter','latex');
    t.Color = 'white';
    t.FontSize = 25;
    t.FontWeight = 'bold';
    t.EdgeColor = 'white';
    t.FitBoxToText = 'on';


    % Print the legend last so the algorithm doesn't place it in front
    % of the text box
    legend('$F_{\uparrow}/F_0$ Monte Carlo','$F_{\downarrow}/F_0$ Monte Carlo',...
        'Interpreter','latex','Location','best')



    s2 = subplot(1,2,2);
    plot(inputs.ssa,linspace(tau_lower_limit, tau_upper_limit, inputs.N_layers))
    hold on
    plot(inputs.g,linspace(tau_lower_limit, tau_upper_limit, inputs.N_layers))
    grid on; grid minor
    set(gca,'YDir','reverse')
    ylabel('$\tau$','Interpreter','latex');
    title({'Multi-layer properties'}, 'Interpreter','latex')
    legend('$\tilde{\omega}$','$g$',...
        'Interpreter','latex','Location','best')

    % Now shift the positions of the subfigures
    s1.Position = s1.Position.*[1 1 1.5 1];

    s2.Position = s2.Position.*[1.25 1 0.5 1];


end



end