%% Plot reflectance curves of constant droplet radius


% By Andrew J. Buggee
%%

function f = plotReflectanceCurves_singleBand(R,inputs,pixels2use, modis)

bands2run = inputs.bands2run;
bands2plot = inputs.bands2plot;
re = inputs.re;
tau_c = inputs.tau_c;
num_pixels = inputs.pixels.num_2calculate;

% extract pixel geometry
sza = pixels2use.res1km.geometry.sza; % solar zenith angle
saz = pixels2use.res1km.geometry.saz; % solar azimuth angle
vza = acosd(pixels2use.res1km.geometry.umu); % viewing zenith angle
vaz = pixels2use.res1km.geometry.phi; % viewing azimuth angle



if num_pixels <=3

    for pp = 1:num_pixels

        f(pp) = figure;

        for bp = 1:size(bands2plot,1)



            for bb = 1:size(bands2plot,2)

                % find indices for bands 2 plot
                bands2plot_index = bands2plot(bp,bb) == bands2run;

                subplot(1,size(bands2plot,2),bb);

                bandVals = modisBands(bands2plot(bp,bb));

                lgnd_str = cell(1,length(re));

                for rr = 1:length(re)


                    reflectance = R(pp,rr,:,bands2plot_index);
                    plot(tau_c,reflectance(:));
                    hold on
                    lgnd_str{rr} = ['r_{e} = ',num2str(re(rr)),' \mum'];
                    set(gcf, 'Position', [0 0 1500 600])

                end

                % plot line of constant reflectance for the MODIS
                % measurement for each band plotted
                [r,c] = ind2sub(pixels2use.res1km.size, pixels2use.res1km.index(pp));
                yline(modis.EV1km.reflectance(r, c, bands2plot(bb)), ...
                    'LineWidth', 4, 'Linestyle','--', 'Color', 'black','Label','MODIS Measurement',...
                    'FontSize', 20, 'Interpreter','latex')


                title([num2str(bandVals(bp,bb)),' nm'])
                xlabel('\tau_{c}'); ylabel('Reflectance')
                grid on; grid minor

                if bb == 1
                    legend(lgnd_str,'Location','best')

                    dim = [0.477900552486189 0.041953385127636 0.148543445504772 0.0305216426193119];
                    str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
                        ' vaz = ',num2str(vaz(pp))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',14, 'EdgeColor', 'w');
                end


            end

        end


    end


elseif num_pixels > 3

    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    rand_index = randsample(num_pixels,3); % random sampling without replacement

    % define set of colors
    C = mySavedColors(1:length(re), 'fixed');

    for pp = 1:length(rand_index)

        for bp = 1:size(bands2plot,1)

            f(pp) = figure;

            for bb = 1:size(bands2plot,2)

                % find indices for bands 2 plot
                bands2plot_index = bands2plot(bp,bb) == bands2run;

                subplot(1,size(bands2plot,2),bb);

                bandVals = modisBands(bands2plot(bp,bb));

                lgnd_str = cell(1,length(re));

                for rr = 1:length(re)

                    reflectance = R(rand_index(pp),rr,:,bands2plot_index);
                    plot(tau_c,reflectance(:), 'color',C(rr,:));
                    hold on
                    lgnd_str{rr} = ['$r_{e} = $',num2str(re(rr)),' $\mu m$'];
                    set(gcf, 'Position', [0 0 1500 600])


                end

                % plot line of constant reflectance for the MODIS
                % measurement for each band plotted
                [r,c] = ind2sub(pixels2use.res1km.size, pixels2use.res1km.index(rand_index(pp)));
                yline(modis.EV1km.reflectance(r, c, bands2plot(bb)), ...
                    'LineWidth', 4, 'Linestyle','--', 'Color', 'black','Label','MODIS Measurement',...
                    'FontSize', 20, 'Interpreter','latex')


                if bb==1
                    % plot line of constant optical depth representing the
                    % MODIs retrieved value
                    [r,c] = ind2sub(pixels2use.res1km.size, pixels2use.res1km.index(rand_index(pp)));
                    xline(modis.cloud.optThickness17(r, c, bands2plot(bb)), ...
                        'LineWidth', 4, 'Linestyle','--', 'Color', 'black','Label','MODIS $\tau_c$',...
                        'FontSize', 20, 'Interpreter','latex')

                    % Print legend
                    legend(lgnd_str,'Location','best', 'Interpreter', 'latex', 'Fontsize',20)
                    
                    % print text box defining pixel index and pixel
                    % geometry
                    dim = [0.35 0.03 0.148543445504772 0.0305216426193119];
                    str = ['Pixel-idx = ', num2str(pixels2use.res1km.index(rand_index(pp))),...
                        '  sza = ',num2str(sza(rand_index(pp))),' saz = ',num2str(saz(rand_index(pp))),...
                        ' vza = ',num2str(vza(rand_index(pp))),' vaz = ',num2str(vaz(rand_index(pp)))];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',16, 'EdgeColor','w', 'Interpreter','latex');

                elseif bb==2
                    % plot line representing the MODIS retrieved droplet
                    % size
                    [r,c] = ind2sub(pixels2use.res1km.size, pixels2use.res1km.index(rand_index(pp)));
                    reflectance = reshape(R(rand_index(pp),:,:,bands2plot_index), length(re), length(tau_c));
                    [tau_mat, re_mat] = meshgrid(tau_c, re);
                    modis_re_interp = interp2(tau_mat, re_mat, reflectance,...
                        tau_c, linspace(modis.cloud.effRadius17(r,c),modis.cloud.effRadius17(r,c),length(tau_c)));

                    plot(tau_c,modis_re_interp,'LineWidth', 4, 'Linestyle','--', 'Color', 'black')

                    % Print legend for the new interpolated MODIS re curve
                    % skip the legend entries that were used in the other
                    % panel and account for the MODIS measurement line
                    new_str = repmat("",1,length(re)+1);
                    % add the new string legend
                    new_str(end+1) = 'MODIS $r_e$';
                    legend(new_str,'Location','best', Interpreter='latex')


                end



                title(['LibRadTran modeled reflectance at ',num2str(bandVals(1)),' nm'], 'Interpreter','latex')
                xlabel('$\tau_{c}$', 'Interpreter', 'latex','FontSize', 24)
                ylabel('Reflectance', 'Interpreter', 'latex','FontSize', 24)
                grid on; grid minor




            end

        end

    end





end




