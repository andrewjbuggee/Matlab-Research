%% ----- Plot two bands of the Reflectance Function against one another -----

% ----- INPUTS -----

% kingORsurf - a string entry that tells the function to create one of two
% plots:
%       (1) the traditional 2D king and nakajima plot where the reflectance
%       can be seen for two bands to vary with droplet size and optical
%       depth
%       (2) a 3D surface plot which shows how small droplets cause
%       redundancy of reflectance, resulting in more than 1 possible
%       solution

% By Andrew J. Buggee

%%

function [] = nakajima_king_reflectance_plot_HySICS(simulated_measurement, modeled_reflectance, inputs, kingORsurf)

% extract inputs

re = inputs.RT.re;
tau_c = inputs.RT.tau_c;

bands2run = inputs.bands2run;
bands2plot = inputs.bands2plot;



if length(bands2plot)~=2

    error('Can only plot two bands at a time')

end


band_vals = mean(inputs.RT.wavelengths2run, 2);



% reshape the libRadTran forward modeled reflectance to be an array where one
% dimension varies with optical depth and another varies with
% effective radius
modelRefl_band1 = reshape(modeled_reflectance(1:2:length(modeled_reflectance)),...
    length(tau_c), length(re));
modelRefl_band2 = reshape(modeled_reflectance(2:2:length(modeled_reflectance)),...
    length(tau_c), length(re));

% Create a legend string
legend_str = cell(1, length(tau_c) + 1 + length(re));
% set the first length(tau_c)+1 entries to be empty strings
for ss = 1:length(tau_c)+1
    legend_str{ss} = '';
end


% R is a 4-D matrix.
%   dim(R,1) = changing pixels
%   dim(R,2) = changing effective radius
%   dim(R,3) = changing optical thickness
%   dim(R,4) = changing spectral band



% Plot values with constant particle radius, and vary the optical thickness



if strcmp(kingORsurf, 'king')==true



        figure;

        % Step through each effective radius. Each line represents the reflectance
        % at a constant droplet size, with varrying optical thickness.
        for rr = 1:length(re)

            plot(modelRefl_band1(:,rr), modelRefl_band2(:,rr),...
                '.-', 'MarkerSize',50,'LineWidth',1.5);
            hold on

            % store legend string for later
            legend_str{length(tau_c) + rr} = ['$r_e = $', num2str(re(rr))];

        end

        % set up color order for each curve
        colororder(mySavedColors(1:length(re),'fixed'));

        % Now plot the simulated HySICS measurement on top
        plot(simulated_measurement.Refl_model(inputs.bands2run_from_set_of_measurements(1)),...
            simulated_measurement.Refl_model(inputs.bands2run_from_set_of_measurements(2)), 'x',...
            'MarkerSize',10, 'Color','black');


        % ------ Plot lines of constant optical depth ------
        % Now step through each optical depth. Each line represents the reflectance
        % at a constant optical depth, with varrying effective radius.
        for tt = 1:length(tau_c)

            x = modelRefl_band1(tt,:);
            y = modelRefl_band2(tt,:);

            t = plot(x, y, 'LineStyle','--', 'Color','k');

            % add line label on plot
            if tt==1
                text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",25, "FontWeight","bold")
                hold on
            else
                text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",25, "FontWeight","bold")
                hold on
            end

            % place the dotted line below the lines of constant radius
            uistack(t, 'bottom');

        end

        % Add text to indicate the black lines are lines of constant optical
        % thickness
        text(1.05*x(end), 0.85*y(end), '$\tau_c$','Interpreter','latex',"FontSize",25, "FontWeight","bold")



        % set up plot stuff
        grid on; grid minor
        xlabel(['Reflectance ', num2str(round(band_vals(1))), ' $nm$'],Interpreter='latex')
        ylabel(['Reflectance ', num2str(round(band_vals(2))), ' $nm$'],Interpreter='latex')

        % Create textbox
    annotation('textbox',[0.3 0.8 0.131464712269273 0.0623003194888175],...
        'String',{['vza = ',num2str(simulated_measurement.inputs.RT.vza), '$^{\circ}$  sza = ',...
        num2str(simulated_measurement.inputs.RT.sza), '$$^{\circ}$']}, ...
        'FontSize', 20,...
        'Interpreter','latex',...
        'FitBoxToText','on');



        % set the last string entry to be MODIS value
        legend_str{end} = 'HySICS';

        legend(legend_str, 'Interpreter','latex','Location','best' , 'FontSize', 20, 'FontWeight','bold')
        title('Simulated HySICS Reflectance for a vertically homogeneous cloud','Interpreter','latex')
        set(gcf,"Position", [0 0 1300 800])









elseif strcmp(kingORsurf, 'surf')==true

    re_mat = repmat(re', 1, length(tau_c));
    tau_c_mat = repmat(tau_c, length(re), 1);

    % Plot across pixels
    for pp = 1:num_pixels_2run


        % plot the visible band first
        f =figure;
        subplot(1,2,1)
        % plot the first band
        surf(tau_c_mat, re_mat, reshape(modeled_reflectance(index_2plot(pp),:,:,index_bands(1)),length(re),[]));
        
        % interpolate between points to smooth the surface
        shading interp

        % set up plot stuff
        grid on; grid minor
        ylabel('$r_e$ ($\mu m$)', Interpreter='latex')
        xlabel('$\tau_c$' ,Interpreter='latex')


        title(['Reflectance ', num2str(band_vals(1,1)), ' $nm$'],Interpreter='latex')



        % Create textbox
        annotation('textbox',[0.4 0.8 0.131464712269273 0.0623003194888175],...
            'String',{['vza = ',num2str(simulated_measurement.sensor.zenith(index_2plot(pp))), '$^{\circ}$  sza = ',...
            num2str(simulated_measurement.solar.zenith(index_2plot(pp))), '$$^{\circ}$', ...
            'pixel idx = ', num2str(pixels2use.res1km.index(index_2plot(pp)))]},...
            'FontSize', 20,...
            'Interpreter','latex',...
            'FitBoxToText','on');


        subplot(1,2,2)
        % plot the second band
        surf(tau_c_mat, re_mat, reshape(modeled_reflectance(index_2plot(pp),:,:,index_bands(2)),length(re),[]));

        % interpolate between points to smooth the surface
        shading interp

        % set up plot stuff
        grid on; grid minor
        ylabel('$r_e$ ($\mu m$)', Interpreter='latex')
        xlabel('$\tau_c$' ,Interpreter='latex')


        title(['Reflectance ', num2str(band_vals(2,1)), ' $nm$'],Interpreter='latex')


        set(gcf,"Position", [0 0 1300 800])



    end


else

    error([newline, 'I dont recognize the kingORsurf string input', newline,...
        'Only two plot options: the King and Najkajima 2D plot or a surface plot',newline])

end





end


