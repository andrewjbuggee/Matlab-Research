%% ----- Plot two bands of the Reflectance Function against one another -----



% By Andrew J. Buggee

%%

function [] = plot_estimated_reflectance_with_modis(modis,R,modisInputs, pixels2use)

% extract inputs

re = modisInputs.re;
tau_c = modisInputs.tau_c;
pixel_row = pixels2use.res1km.row;
pixel_col = pixels2use.res1km.col;

modis_band_vals = modisBands([modisInputs.bands2run]);


% Figure out the number of subfigures there are, and the number of rows and
% columns ther should be
num_subFigures = length(modisInputs.bands2run);


if num_subFigures<=4

    % In this case, we will only have 1 row
    num_figure_cols = num_subFigures;
    num_figure_rows = 1;

elseif rem(num_subFigures, 2)==1

    % In this case, we will have 2 rows
    num_figure_cols = (num_subFigures+1)/2;
    num_figure_rows = 2;

end



% Plot values with constant optical thickness, and varry the effective
% droplet size

% Create legend string
legend_str =[string(tau_c), 'MODIS'];
% set up the first legend entry with the proper variable
for ii = 1:length(tau_c)
    legend_str(ii) = ['$\tau_c = $', legend_str{ii}];
end

% R is a 4-D matrix.
%   dim(R,1) = changing pixels
%   dim(R,2) = changing effective radius
%   dim(R,3) = changing optical thickness
%   dim(R,4) = changing spectral band

% create a new figure for each pixel

for pp = 1:size(R,1)

    figure;

    % Now step through each subpanel, which is a different wavelength
    for bb = 1:num_subFigures

        subplot(num_figure_rows, num_figure_cols, bb)

        plot(re, reshape(R(pp, :, :, bb),length(re), length(tau_c)), '.-', 'MarkerSize',15,'LineWidth',1);
        grid on; grid minor
        xlabel('$r_e$  ($\mu m$)','Interpreter','latex')
        title([num2str(modis_band_vals(bb,1)), ' $nm$'], 'Interpreter','latex', 'FontSize',20)

        
        hold on;
        % Now plot the MODIS reflectance measurements as horizontal lines
        yline(modis.EV1km.reflectance(pixel_row(pp), pixel_col(pp),modisInputs.bands2run(bb)), '--',...
            'color','black','LineWidth',1)


        % not everything needs to be printed in each panel
        if bb==1
            ylabel('Reflectance ($1/sr$)','Interpreter','latex')

        elseif bb==3
            legend(legend_str,'Interpreter','latex', 'Location','best')
        end

        % Define the color order
        colororder(mySavedColors(1:length(tau_c), 'fixed'))

        hold on
    end

    % Include a super title
    sgtitle('Comparing Model Reflectances with MODIS','interpreter','latex', 'Fontsize',30)
    % set size and position
    set(gcf, 'Position', [0 0 1250 500])


end












end


