%% ----- Plot Earth Viewing Spectral Data from MODIS L1B data -----

% inputs are: 
% (data,sensor_parameters,location_params,plot_resolution,plot_quantity)

% By Andrew J. Buggee
%%

function [] = plotL1B_data(EV,geo,plotRes,quantity)


% ---- define global parameters ----

nbins = 100;  % number of deiscrete bins for histogram data
pauseTime = 0.33; % seconds

% geolocation values are at 1km resolution, always
xBoundries = [ceil(geo.long(end,1)),floor(geo.long(end,end))]; % Boundries of the image in terms of longitude
yBoundries = [floor(geo.lat(1,1)),ceil(geo.lat(end,1))]; % boundaries of the image in terms of latitude



% ----- create plots -----


if strcmp(plotRes,'250')
    
    % interpolate lat and long so they are the same dimensions as the science
    % data
    
    geoX = geo.long(end,:);
    x1 = 1:length(geoX); % create a simple independent variable
    x2 = linspace(1,length(geoX),size(EV.m250.radiance,2)); % new independent variable with the same spread but finer values
    geoX = interp1(x1,geoX,x2);
    
    geoY = geo.lat(:,1);
    y1 = 1:length(geoY); % create a simple independent variable
    y2 = linspace(1,length(geoY),size(EV.m250.radiance,1)); % new independent variable with the same spread but finer values
    geoY = interp1(y1,geoY,y2);
    
    
    
    if strcmp(quantity,'radiance')
        
        
        for ii = 1:size(EV.m250.radiance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX,fliplr(geoY),EV.m250.radiance(:,:,ii));
            colorbar
            title('Radiance (W/m^{2}/\mum/sr)')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m250.radiance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m250.bands.center(ii)),' \mum'])
            xlabel('Radiance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    elseif strcmp(quantity,'reflectance')
        
        for ii = 1:size(EV.m250.reflectance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX,fliplr(geoY),EV.m250.reflectance(:,:,ii));
            colorbar
            title('Reflectance')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m250.reflectance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m250.bands.center(ii)),' \mum'])
            xlabel('Reflectance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    else
        error('Not a valid quantity to plot')
        
    end
    
    
elseif strcmp(plotRes,'500')
    
        % interpolate lat and long so they are the same dimensions as the science
    % data
    
    geoX = geo.long(end,:);
    x1 = 1:length(geoX); % create a simple independent variable
    x2 = linspace(1,length(geoX),size(EV.m500.radiance,2)); % new independent variable with the same spread but finer values
    geoX = interp1(x1,geoX,x2);
    
    geoY = geo.lat(:,1);
    y1 = 1:length(geoY); % create a simple independent variable
    y2 = linspace(1,length(geoY),size(EV.m500.radiance,1)); % new independent variable with the same spread but finer values
    geoY = interp1(y1,geoY,y2);
    
    
    if strcmp(quantity,'radiance')
        
        
        for ii = 1:size(EV.m500.radiance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX,fliplr(geoY),EV.m500.radiance(:,:,ii));
            colorbar
            title('Radiance (W/m^{2}/\mum/sr)')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m500.radiance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m500.bands.center(ii)),' \mum'])
            xlabel('Radiance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    elseif strcmp(quantity,'reflectance')
        
        for ii = 1:size(EV.m500.reflectance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX,fliplr(geoY),EV.m500.reflectance(:,:,ii));
            colorbar
            title('Reflectance')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m500.reflectance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m500.bands.center(ii)),' \mum'])
            xlabel('Reflectance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    else
        error('Not a valid quantity to plot')
        
    end
    
elseif strcmp(plotRes,'all')
    
        
    geoX = geo.long(end,:);
    x1 = 1:length(geoX); % create a simple independent variable
    x2 = linspace(1,length(geoX),size(EV.m250.radiance,2)); % new independent variable with the same spread but finer values
    geoX_250 = interp1(x1,geoX,x2);
    x2 = linspace(1,length(geoX),size(EV.m500.radiance,2)); % new independent variable with the same spread but finer values
    geoX_500 = interp1(x1,geoX,x2);
    
    
    geoY = geo.lat(:,1);
    y1 = 1:length(geoY); % create a simple independent variable
    y2 = linspace(1,length(geoY),size(EV.m250.radiance,1)); % new independent variable with the same spread but finer values
    geoY_250 = interp1(y1,geoY,y2);
    y2 = linspace(1,length(geoY),size(EV.m500.radiance,1)); % new independent variable with the same spread but finer values
    geoY_500 = interp1(y1,geoY,y2);
        
    
    if strcmp(quantity,'radiance')
        
        
        for ii = 1:size(EV.m250.radiance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX_250,fliplr(geoY_250),EV.m250.radiance(:,:,ii));
            colorbar
            title('Radiance (W/m^{2}/\mum/sr)')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m250.radiance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m250.bands.center(ii)),' \mum'])
            xlabel('Radiance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
        
        for ii = 1:size(EV.m500.radiance,3)
            
            figure;
            subplot(1,2,1)
            imagesc(geoX_500,fliplr(geoY_500),EV.m500.radiance(:,:,ii));
            colorbar
            title('Radiance (W/m^{2}/\mum/sr)')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m500.radiance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m500.bands.center(ii)),' \mum'])
            xlabel('Radiance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    elseif strcmp(quantity,'reflectance')
        
        for ii = 1:size(EV.m250.reflectance,3)
            figure;
            subplot(1,2,1)
            imagesc(geoX_250,fliplr(geoY_250),EV.m250.reflectance(:,:,ii));
            colorbar
            title('Reflectance')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m250.reflectance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m250.bands.center(ii)),' \mum'])
            xlabel('Reflectance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
        
        for ii = 1:size(EV.m500.reflectance,3)
            
            figure;
            subplot(1,2,1)
            imagesc(geoX_500,fliplr(geoY_500),EV.m500.reflectance(:,:,ii));
            colorbar
            title('Reflectance')
            ylabel('Latitude')
            xlabel('Longitude')
            
            xticks([xBoundries(1),xBoundries(2)]);
            xticklabels({num2str(xBoundries(1)),num2str(xBoundries(2))});
            
            yticks([yBoundries(2),yBoundries(1)]);
            yticklabels({num2str(yBoundries(1)),num2str(yBoundries(2))});
            pause(pauseTime)
            
            % plot the histogram of counts
            
            subplot(1,2,2)
            histogram(EV.m500.reflectance(:,:,ii),nbins)
            title(['Band: ',num2str(EV.m500.bands.center(ii)),' \mum'])
            xlabel('Reflectance')
            ylabel('Frequency')
            set(gca,'yscale','log')
        end
        
    end
    
    
else
    
    error('Not a valid plotFlag')
    
end



end
