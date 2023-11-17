%% ----- Plot two bands of the Reflectance Function against one another -----

% ----- INPUTS -----

% kingORsurf - a string entry that tells the function to create one of two
% plots:
%       (1) the traditional 2D king an nakajima plot where the reflectance
%       can be seen for two bands to vary with droplet size and optical
%       depth
%       (2) a 3D surface plot which shows how small droplets cause
%       redundancy of reflectance, resulting in more than 1 possible
%       solution

% By Andrew J. Buggee

%%

function [] = plot2ReflectanceFuncBands(modis,R,modisInputs, pixels2use, kingORsurf)

% extract inputs

re = modisInputs.RT.re;
tau_c = modisInputs.RT.tau_c;
pixel_row = pixels2use.res1km.row;
pixel_col = pixels2use.res1km.col;
bands2run = modisInputs.bands2run;
bands2plot = modisInputs.bands2plot;



if length(bands2plot)~=2

    error('Can only plot two bands at a time')

end


modis_band_vals = modisBands([modisInputs.bands2plot]);

% find the indices needed to plot the bands listed above
for ii = 1:length(bands2plot)
    index_bands(ii) = find(bands2run==bands2plot(ii));
end


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

% Only plot 3 pixels at a time
if size(R,1)>3

    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    index_2plot = randsample(size(R,1),3); % random sampling without replacement

    % 3 examples will be run
    num_pixels_2run = 3;

else

    index_2plot = 1:size(R,1);
    num_pixels_2run = size(R,1);

end


if strcmp(kingORsurf, 'king')==true

    % Plot across pixels
    for pp = 1:num_pixels_2run

        figure;

        % Step through each effective radius. Each line represents the reflectance
        % at a constant droplet size, with varrying optical thickness.
        for rr = 1:length(re)

            plot(reshape(R(index_2plot(pp),rr,:,index_bands(1)),1,[]), reshape(R(index_2plot(pp),rr,:,index_bands(2)),1,[]),...
                '.-', 'MarkerSize',50,'LineWidth',1.5);
            hold on

            % store legend string for later
            legend_str{length(tau_c) + rr} = ['$r_e = $', num2str(re(rr))];

        end

        % set up color order for each curve
        colororder(mySavedColors(1:length(re),'fixed'));

        % Now plot the MODIS measurement on top
        plot(modis.EV1km.reflectance(pixel_row(index_2plot(pp)), pixel_col(index_2plot(pp)),modisInputs.bands2run(index_bands(1))),...
            modis.EV1km.reflectance(pixel_row(index_2plot(pp)), pixel_col(index_2plot(pp)),modisInputs.bands2run(index_bands(2))), 'x',...
            'MarkerSize',10, 'Color','black');


        % ------ Plot lines of constant optical depth ------
        % Now step through each optical depth. Each line represents the reflectance
        % at a constant optical depth, with varrying effective radius.
        for tt = 1:length(tau_c)

            x = reshape(R(index_2plot(pp),:,tt,index_bands(1)),1,[]);
            y = reshape(R(index_2plot(pp),:,tt,index_bands(2)),1,[]);

            t = plot(x, y, 'LineStyle','--', 'Color','k');

            % add line label on plot
            if tt==1
                text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",20, "FontWeight","bold")
                hold on
            else
                text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",20, "FontWeight","bold")
                hold on
            end

            % place the dotted line below the lines of constant radius
            uistack(t, 'bottom');

        end

        % Add text to indicate the black lines are lines of constant optical
        % thickness
        text(1.05*x(end), 0.85*y(end), '$\tau_c$','Interpreter','latex',"FontSize",20, "FontWeight","bold")



        % set up plot stuff
        grid on; grid minor
        xlabel(['Reflectance ', num2str(modis_band_vals(1,1)), ' $nm$'],Interpreter='latex')
        ylabel(['Reflectance ', num2str(modis_band_vals(2,1)), ' $nm$'],Interpreter='latex')

        % Create textbox
    annotation('textbox',[0.3 0.8 0.131464712269273 0.0623003194888175],...
        'String',{['vza = ',num2str(modis.sensor.zenith(index_2plot(pp))), '$^{\circ}$  sza = ',...
        num2str(modis.solar.zenith(index_2plot(pp))), '$$^{\circ}$', ...
        '  pixel-idx = ', num2str(pixels2use.res1km.index(index_2plot(pp)))]},...
        'FontSize', 20,...
        'Interpreter','latex',...
        'FitBoxToText','on');



        % set the last string entry to be MODIS value
        legend_str{end} = 'MODIS';

        legend(legend_str, 'Interpreter','latex','Location','best' , 'FontSize', 20, 'FontWeight','bold')
        title('Simulated Reflectance','Interpreter','latex')
        set(gcf,"Position", [0 0 1300 800])





    end





elseif strcmp(kingORsurf, 'surf')==true

    re_mat = repmat(re', 1, length(tau_c));
    tau_c_mat = repmat(tau_c, length(re), 1);

    % Plot across pixels
    for pp = 1:num_pixels_2run


        % plot the visible band first
        f =figure;
        subplot(1,2,1)
        % plot the first band
        surf(tau_c_mat, re_mat, reshape(R(index_2plot(pp),:,:,index_bands(1)),length(re),[]));
        
        % interpolate between points to smooth the surface
        shading interp

        % set up plot stuff
        grid on; grid minor
        ylabel('$r_e$ ($\mu m$)', Interpreter='latex')
        xlabel('$\tau_c$' ,Interpreter='latex')


        title(['Reflectance ', num2str(modis_band_vals(1,1)), ' $nm$'],Interpreter='latex')



        % Create textbox
        annotation('textbox',[0.4 0.8 0.131464712269273 0.0623003194888175],...
            'String',{['vza = ',num2str(modis.sensor.zenith(index_2plot(pp))), '$^{\circ}$  sza = ',...
            num2str(modis.solar.zenith(index_2plot(pp))), '$$^{\circ}$', ...
            'pixel idx = ', num2str(pixels2use.res1km.index(index_2plot(pp)))]},...
            'FontSize', 20,...
            'Interpreter','latex',...
            'FitBoxToText','on');


        subplot(1,2,2)
        % plot the second band
        surf(tau_c_mat, re_mat, reshape(R(index_2plot(pp),:,:,index_bands(2)),length(re),[]));

        % interpolate between points to smooth the surface
        shading interp

        % set up plot stuff
        grid on; grid minor
        ylabel('$r_e$ ($\mu m$)', Interpreter='latex')
        xlabel('$\tau_c$' ,Interpreter='latex')


        title(['Reflectance ', num2str(modis_band_vals(2,1)), ' $nm$'],Interpreter='latex')


        set(gcf,"Position", [0 0 1300 800])



    end


else

    error([newline, 'I dont recognize the kingORsurf string input', newline,...
        'Only two plot options: the King and Najkajima 2D plot or a surface plot',newline])

end





end


