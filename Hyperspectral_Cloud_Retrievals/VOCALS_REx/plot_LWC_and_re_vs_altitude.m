%% Plot the Liquid Water Content and the Effective Radius versus Altitude 


% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_vs_altitude(vert_profiles, indices)


% Plot the number of curves within the vert_profiles structure

N_cuvres = length(indices);

figure;
for nn = 1:N_cuvres
    
    % First plot the LWC
    subplot(1,2,1); plot(vert_profiles.lwc{indices(nn)}, vert_profiles.altitude{indices(nn)})
    hold on

    % next plot the effective radius
    subplot(1,2,2); plot(vert_profiles.re{indices(nn)}, vert_profiles.altitude{indices(nn)}); 
    hold on


end

% Make each subplot pretty
subplot(1,2,1)
grid on; grid minor; 
xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
ylabel('Altitude ($m$)', 'Interpreter','latex');


subplot(1,2,2)
grid on; grid minor; 
xlabel('re ($\mu m$)', 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1000 550])




end