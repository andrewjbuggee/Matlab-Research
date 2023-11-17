%% Plot all scattering events within our medium


function plotScatteringEvents(photon_tau_position, tau_y_upperLimit)

figure;
xline(0,'LineWidth',0.5,'Color','black')
hold on;
yline(0,'LineWidth',1.5,'Color','black')
hold on
yline(tau_y_upperLimit,'LineWidth',1.5,'Color','black')
hold on;
grid on; grid minor

% Set figure position and size


ylim([-1,tau_y_upperLimit+1])
xlim([min(photon_tau_position(:,1))-1,max(photon_tau_position(:,1))+1])
title('Photon position','Interpreter','latex')
xlabel('$\tau_x$','Interpreter','latex')
ylabel('$\tau_y$','Interpreter','latex')

set(gca,'YDir','reverse')

drawnow



for nn = 1:size(photon_tau_position,1)-1


    plot([photon_tau_position(nn,1), photon_tau_position(nn+1,1)],...
              [photon_tau_position(nn,2), photon_tau_position(nn+1,2)]);
    hold on
    drawnow
    pause(0.25)




end







end


