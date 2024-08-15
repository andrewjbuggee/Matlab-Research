% By Andrew John Buggee

clear variables
%% ---- TWO STREAM MODEL ---

% This code assumes the following:
%   (1) a semi-infinite cloud, with a hard boundary at the top
%   and extending to infinity in the other direction.
%   (2) We assume photons can only travel up or down. Therefore the two end
%   states are (a) the photon is absorbed in the cloud or (b) the photon
%   scatters back out the cloud top

% Bohren and Clothiaux list a solution to this infinite two-stream problem
% in their book, starting on page 264. The reflection out the cloud top,
% the irradiance that we measure, is defined on page 265, and called R.


% Compute the reflectance out of the cloud top across many wavelengths. We
% need the asymmetry parameter and the absorption coefficient

% compute the single scattering albedo and asymmetry parameter using mie
% theory for some wavelength w and some radius r. Also, specify whether to
% use a monodispersed distribution or a gamma droplet distribution

droplet_distribution = 'gamma';
effective_radius = 10;                        % microns - effective droplet radii

% --- wavelength limits are between [100,3000] nanometers ---
%wavelength = 100:10:3000;                      % nm - wavelength
wavelength = 350:100:2300;                       % nm - wavelength range of the HySICS spectrometer
%wavelength = [1000, 1100, 1200];                               % nm - wavelength

mie_properties = zeros(length(wavelength),8,length(effective_radius));
legend_str = cell(1,length(effective_radius));


for rr = 1:length(effective_radius)
    mie_properties(:,:,rr) = interp_mie_computed_tables([wavelength',repmat(effective_radius(rr),length(wavelength),1)],droplet_distribution,false);

    legend_str{rr} = ['$r_{e}$ = ',num2str(effective_radius(rr)),' $\mu m$'];
end

single_scattering_albedo = reshape(mie_properties(:,6,:),length(wavelength),[]);
g = reshape(mie_properties(:,7,:), length(wavelength),[]);

%% Compute the Two Stream Reflectivity


%complex_index_refraction = mie_properties(:,4,:);

%absoprtion_coefficient = 4*pi*complex_index_refraction./(wavelength'./1e3);         % microns^(-1)

% ----- using the Towmey and Bohren approximation for absorption ---
% R = (sqrt(1 - g.*(1 - 0.85*effective_radius*absoprtion_coefficient)) - sqrt(0.85*effective_radius*absoprtion_coefficient))./...
%     (sqrt(1 - g.*(1 - 0.85*effective_radius*absoprtion_coefficient)) + sqrt(0.85*effective_radius*absoprtion_coefficient));


% ----- Using the mie estimated values for single scattering ----
% R_top is the amount of irradiance reflected out the top of the cloud
R_top = (sqrt(1 - g.*single_scattering_albedo) - sqrt(1 - single_scattering_albedo))./...
    (sqrt(1 - g.*single_scattering_albedo) + sqrt(1 - single_scattering_albedo));



% plot the results of reflectivity versus wavelength for several
% different droplet sizes

if length(wavelength)>length(effective_radius)

    figure; plot(wavelength, R_top)
    xlabel('Wavelength (nm)','Interpreter','latex')
    ylabel('Reflectivity','Interpreter','latex')
    title(['Two-Stream Reflectivity: semi-infinite ', droplet_distribution,'-dispersed cloud'],'Interpreter','latex')
    grid on; grid minor
    legend(legend_str, 'Interpreter','latex','Location','best','FontSize',18)
    set(gcf,'Position',[0 0 1200, 800])

elseif length(wavelength)==1 && length(effective_radius)>1

    figure; plot(effective_radius, R_top)
    xlabel('Effective Radius $(\mu m)$','Interpreter','latex')
    ylabel('Reflectivity','Interpreter','latex')
    title(['Two-Stream Reflectivity: semi-infinite ', droplet_distribution,'-dispersed cloud'],'Interpreter','latex')
    grid on; grid minor
    legend([num2str(wavelength),' nm'], 'Interpreter','latex','Location','best','FontSize',18)
    set(gcf,'Position',[0 0 1200, 800])

end



%% Compare the 2-stream reflectivity with the approximations by Twomey and Bohren


% plot the results of reflectivity versus wavelength for several
% different droplet sizes

% plot colors
randomColors = rand(4, 3);

solar_zenithAngle = 0;


if length(wavelength)>length(effective_radius)

    R_TB = 1 - sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo)) .* H(cosd(solar_zenithAngle), single_scattering_albedo,0)';


    % Lets aslo compute the simple linear approximation for absorption by a
    % single water droplet in the geometric optics limit
    absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)

    R_TB_simple = 1 - H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* sqrt((0.85 * absorption_coefficient .* effective_radius)./(1 - g.*(1 - 0.85 * effective_radius.*absorption_coefficient)));

    R_TB_superSimple = 1 - 2.2 * H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* sqrt(absorption_coefficient .* effective_radius);



    figure; plot(wavelength, R_top,'-','Color',randomColors(1,:))
    hold on;
    plot(wavelength, R_TB,'--','Color',randomColors(2,:))
    plot(wavelength, R_TB_simple,'-.','Color',randomColors(3,:))
    plot(wavelength, R_TB_superSimple, 'LineStyle',':','Color',randomColors(4,:))
    xlabel('Wavelength (nm)','Interpreter','latex')
    ylabel('Reflectivity','Interpreter','latex')
    title(['Comparing Reflectivities for a semi-infinite ', droplet_distribution,'-dispersed cloud'],'Interpreter','latex')
    grid on; grid minor
    legend('2-Stream','Absorption w/ H-functions','$1 - \tilde{\omega} \approx 0.85 r k$','$\bar{g}$ \& $\bar{k}$','Location','best',...
        'Interpreter','latex')
    set(gcf,'Position',[0 0 1200, 800])
    ylim([-1,1])

elseif length(wavelength)==1 && length(effective_radius)>1

    figure; plot(effective_radius, R_top)
    xlabel('Effective Radius $(\mu m)$','Interpreter','latex')
    ylabel('Reflectivity','Interpreter','latex')
    title(['Two-Stream Reflectivity: semi-infinite ', droplet_distribution,'-dispersed cloud'],'Interpreter','latex')
    grid on; grid minor
    legend([num2str(wavelength),' nm'], 'Interpreter','latex','Location','best','FontSize',18)
    set(gcf,'Position',[0 0 1200, 800])

end


%% Least Squares Solution - Can I retrieve r?

% no, this didn't work. Why not?
% Something about the assumptions of the two equations I'm equalizing must
% be wrong. Check the assumptions of the Two-stream reflectivity I'm using,
% and compare them with those of the Twomey Bohren approximation for
% absorption.

index_wavelength = [2, 50, 80];         % wavelengths to use

absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)


G = [ones(length(index_wavelength),1), -0.85* absorption_coefficient(index_wavelength)];         % linear model

d = (G' * G)^(-1) *G' * R_top(index_wavelength);

%% How does the two stream multiple scattering solutions compare to Twomey and Bohren?

% We solved the two-stream multiple scattering reflectivity of a
% semi-infinite cloud. Twomey and Bohren do the same using H-functions to
% solve for the absorption due to normal incident flux. These should be similar.

radius2plot = 10;
index2plot = effective_radius==radius2plot;

solar_zenithAngle = 0;                  % deg

abs_TB = sqrt((1 - single_scattering_albedo(:,index2plot))./(1 - g(:,index2plot).*single_scattering_albedo(:,index2plot))) .* H(cosd(solar_zenithAngle), single_scattering_albedo(:,index2plot),0)';

% In a semi-infinte layer, only absorption and reflection out the cloud top
% exist as the two photon end states.

R_TB = 1 - abs_TB;

% plot the results of two stream reflectivity and the reflectivity
% estiamted using the absorption approximation by Twomey and Bohren

figure; plot(wavelength, R_top(:,index2plot))
hold on;
plot(wavelength, R_TB)
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Reflectivity','Interpreter','latex')
title(['Reflectivity of a semi-infinite cloud w/ $r_e$ = ',num2str(radius2plot),' $\mu m$'],'Interpreter','latex')
grid on; grid minor
legend('Two-Stream','Twomey-Bohren', 'Interpreter','latex','Location','best','FontSize',18)
set(gcf,'Position',[0 0 1200, 800])


%% Plot Twomey and Bohren's absorption approximation for a semi-infinite cloud with respect to wavelength or effective radius

% Let's use the absorption approximation for a semi-infinite cloud layer
% derived by Twomey and Bohren (equation 3, 1980). Is the absorption
% roughly linear with effective radius?


solar_zenithAngle = 0;                  % deg



if length(wavelength)==1 && length(effective_radius)>1


    abs_TB = sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo)) .* H(cosd(solar_zenithAngle), single_scattering_albedo,0);


    % Lets aslo compute the simple linear approximation for absorption by a
    % single water droplet in the geometric optics limit
    absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)

    abs_TB_simple = H(cosd(solar_zenithAngle), single_scattering_albedo,0) .* sqrt((0.85 * absorption_coefficient .* effective_radius)./(1 - g.*(1 - 0.85 * effective_radius.*absorption_coefficient)));

    abs_TB_superSimple = 2.2 * H(cosd(solar_zenithAngle), single_scattering_albedo,0) .* sqrt(absorption_coefficient .* effective_radius);



    % plot the results of two stream reflectivity and the reflectivity
    % estiamted using the absorption approximation by Twomey and Bohren

    figure; plot(effective_radius, abs_TB)
    hold on; plot(effective_radius, abs_TB_simple,'--')
    plot(effective_radius, abs_TB_superSimple, 'LineStyle',':')
    xlabel('Effective Radius $(\mu m)$','Interpreter','latex')
    ylabel('Absorption','Interpreter','latex')
    title(['Absorption by a semi-infinite cloud w/ $\lambda$ = ',num2str(wavelength),' $nm$'],'Interpreter','latex')
    legend('Absorption w/ H-functions','Simple approximation','Super Simple approx','Location','best')
    grid on; grid minor
    set(gcf,'Position',[0 0 1200, 800])


elseif length(wavelength)>1 && length(effective_radius)==1


    abs_TB = sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo)) .* H(cosd(solar_zenithAngle), single_scattering_albedo,0)';


    % Lets aslo compute the simple linear approximation for absorption by a
    % single water droplet in the geometric optics limit
    absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)

    abs_TB_simple = H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* sqrt((0.85 * absorption_coefficient .* effective_radius)./(1 - g.*(1 - 0.85 * effective_radius.*absorption_coefficient)));

    abs_TB_superSimple = 2.2 * H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* sqrt(absorption_coefficient .* effective_radius);


    figure; plot(wavelength, abs_TB)
    hold on; plot(wavelength, abs_TB_simple,'--')
    plot(wavelength, abs_TB_superSimple, 'LineStyle',':')
    xlabel('Wavelength $(nm)$','Interpreter','latex')
    ylabel('Absorption','Interpreter','latex')
    title(['Absorption by a semi-infinite cloud w/ $r_e$ = ',num2str(effective_radius),' $\mu m$ and $\mu_0$ = ',...
        num2str(cosd(solar_zenithAngle))],'Interpreter','latex')
    legend('Absorption w/ H-functions','$1 - \tilde{\omega} \approx 0.85 r k$','$\bar{g}$ \& $\bar{k}$','Location','best')
    grid on; grid minor
    set(gcf,'Position',[0 0 1200, 800])


else

    error('I dont know what to do!')


end


%% Reproducing Twomey and Bohren figure 4 where fractional absorption varies with solar zenith angle

solar_zenithAngle = 0:5:85;                  % deg


% plot colors
randomColors = rand(length(wavelength), 3);


if length(solar_zenithAngle)>1 && length(effective_radius)==1

    abs_TB = repmat(sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo)),1,length(solar_zenithAngle))...
        .* H(cosd(solar_zenithAngle), single_scattering_albedo,0)';


    % Lets aslo compute the simple linear approximation for absorption by a
    % single water droplet in the geometric optics limit
    absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);  % microns^(-1)

    abs_TB_simple = H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* ...
        repmat(sqrt((0.85 * absorption_coefficient .* effective_radius)./(1 - g.*(1 - 0.85 * effective_radius.*absorption_coefficient))), 1, length(solar_zenithAngle));


    abs_TB_superSimple = 2.2 * H(cosd(solar_zenithAngle), single_scattering_albedo,0)' .* ...
        repmat(sqrt(absorption_coefficient .* effective_radius),1,length(solar_zenithAngle));


    figure;
    legend_str = cell(1,length(wavelength));

    for ww = 1:length(wavelength)

        plot(cosd(solar_zenithAngle), abs_TB_superSimple(ww,:), 'LineStyle','-','Color',randomColors(ww,:))
        hold on
        legend_str{ww} = [num2str(wavelength(ww)/1e3), ' $\mu m$'];
    end

    xlabel('$\mu_0$','Interpreter','latex')
    ylabel('Absorption','Interpreter','latex')
    title(['Absorption by a semi-infinite cloud w/ $r_e$ = ',num2str(effective_radius)],'Interpreter','latex')
    legend(legend_str,'Location','best','Interpreter','latex')
    grid on; grid minor
    set(gcf,'Position',[0 0 1200, 800])


else

    error('I dont know what to do!')

end




%% Reproducing Twomey and Bohren Figure 3 - plotting the single scattering absorption approximation


% plot colors
randomColors = rand(length(wavelength), 3);


if length(effective_radius)>1 


    % Lets aslo compute the simple linear approximation for absorption by a
    % single water droplet in the geometric optics limit
    absorption_coefficient = 4*pi*reshape(mie_properties(:,4,:), length(wavelength),[])./(wavelength'./1e3);        % microns^(-1)

    single_scattering_absorption = 0.85 * repmat(effective_radius,length(wavelength),1).*absorption_coefficient;


    figure;
    legend_str = cell(1,length(wavelength));

    for ww = 1:length(wavelength)

        semilogy(effective_radius, single_scattering_absorption(ww,:), 'LineStyle','-','Color',randomColors(ww,:))
        hold on
        legend_str{ww} = [num2str(wavelength(ww)/1e3), ' $\mu m$'];
    end

    xlabel('Effective Radius $(\mu m)$','Interpreter','latex')
    ylabel('$ 1 - \tilde{\omega}$','Interpreter','latex')
    title(['Single Particle Absorption'],'Interpreter','latex')
    legend(legend_str,'Location','best','Interpreter','latex')
    grid on; grid minor
    set(gcf,'Position',[0 0 1200, 800])


else

    error('I dont know what to do!')

end


%% For an isotropic cloud, the fractional absorption only depends on single scattering albedo and the source direction

% If I compute the fractional absorption over this region, does the
% magnitude ever eclipse 1?

solar_zenithAngle = 0;

if length(wavelength)>1


    abs_TB_isotropic = sqrt(1 - single_scattering_albedo) .* H(cosd(solar_zenithAngle), single_scattering_albedo,0)';
    
    abs_TB_anisotropic = sqrt((1 - single_scattering_albedo)./(1 - g.*single_scattering_albedo)) .* H(cosd(solar_zenithAngle), single_scattering_albedo,0)';



    figure; 
    plot(wavelength, abs_TB_isotropic)
    hold on; plot(wavelength, abs_TB_anisotropic)
    xlabel('Wavelength $(nm)$','Interpreter','latex')
    ylabel('Absorption','Interpreter','latex')
    title(['Absorption by a semi-infinite cloud w/ $r_e$ = ',num2str(effective_radius),' $\mu m$ and $\mu_0$ = ',...
        num2str(cosd(solar_zenithAngle))],'Interpreter','latex')
    legend('Isotropic Cloud','Anisotropic Cloud')
    grid on; grid minor
    set(gcf,'Position',[0 0 1200, 800])


end


%% --- Compute the upward irradiance at different tau depths

% K is defined on page 265 of Bohren and Clothiaux
K = sqrt((1 - single_scattering_albedo).*(1 - g.*single_scattering_albedo));

% Compute the upward moving light normalized by the incident light at
% different tau levels
R_up = @(tau) repmat(R_top,1,length(tau)).*exp(-K*tau);

% defind an optical depth vector
tau = 0:0.1:50;                     % optical depth

figure; plot(R_up(tau),tau)
set(gca,'YDir','reverse')
grid on; grid minor
xlabel('$F_{\uparrow}(\tau)/F_{\downarrow}(\tau=0)$','Interpreter','latex')
ylabel('$\tau$','Interpreter','latex')
title('Fraction of incident light moving upwards at different $\tau$','Interpreter','latex')

% let's plot the derivative of the fraction of upward moving light with
% respect to tau
figure; plot(repmat(R_top.*K,1,length(tau)).*exp(-K*tau),tau)
set(gca,'YDir','reverse')
grid on; grid minor
xlabel('$\partial R_{\uparrow}(\tau)/ \partial\tau$','Interpreter','latex')
ylabel('$\tau$','Interpreter','latex')
title('Distribution of upwards moving light w.r.t. $\tau$','Interpreter','latex')
legend(string(wavelength),'location','best','Interpreter','latex')




