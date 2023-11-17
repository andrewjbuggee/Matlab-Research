%% Create Reflectance Contour Plots



% By Andrew J. Buggee
%%

function f = plotReflectanceContours(R,inputs,pixels2use)

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


if num_pixels <= 3
    
    for pp = 1:num_pixels
        
        f(pp) = figure;
        
        for ii = 1:size(bands2plot,1)
            
            
            
            for jj = 1:size(bands2plot,2)
                
                % find indices for bands 2 plot
                bands2plot_index = bands2plot(ii,jj) == bands2run;
                
                bandVals = modisBands(bands2plot(ii,jj));
                subplot(1,size(bands2plot,2),jj)
                
                reflectance = reshape(R(pp,:,:,bands2plot_index),length(re),length(tau_c));
                contourf(tau_c,re,reflectance,'ShowText','on'); colorbar;
                title(['Reflectance - ',num2str(bandVals(1)),' nm'])
                xlabel('\tau_{c}'); ylabel('r_{e} (\mum)')
                set(gcf, 'Position', [0 0 1000 400])

                
            end
            
            
        end
        
        dim = [.5 0 .3 .3];
        str = ['sza = ',num2str(sza(pp)),' saz = ',num2str(saz(pp)),' vza = ',num2str(vza(pp)),...
            ' vaz = ',num2str(vaz(pp))];
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',14, 'EdgeColor', 'w');
        
        
    end
    
    
elseif num_pixels > 3
    
    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    rand_index = randsample(num_pixels,3); % random sampling without replacement
    
    
    for pp = 1:length(rand_index)
        
        f(pp) = figure;
        
        for ii = 1:size(bands2plot,1)
            
            
            
            for jj = 1:size(bands2plot,2)
                
                % find indices for bands 2 plot
                bands2plot_index = bands2plot(ii,jj) == bands2run;
                
                bandVals = modisBands(bands2plot(ii,jj));
                subplot(1,size(bands2plot,2),jj)
                
                reflectance = reshape(R(rand_index(pp),:,:,bands2plot_index),length(re),length(tau_c));
                contourf(tau_c,re,reflectance,'ShowText','on'); colorbar;
                title(['Reflectance - ',num2str(bandVals(1)),' nm'])
                xlabel('\tau_{c}'); ylabel('r_{e} (\mum)')
                set(gcf, 'Position', [0 0 1000 400])

                
            end
            
            
        end
        
        dim = [.5 0 .3 .3];
        str = ['sza = ',num2str(sza(rand_index(pp))),' saz = ',num2str(saz(rand_index(pp))),' vza = ',num2str(vza(rand_index(pp))),...
            ' vaz = ',num2str(vaz(rand_index(pp)))];
        annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k',...
                        'FontWeight','bold','FontSize',14, 'EdgeColor', 'w');
        
        
    end
    
end




end