% ----- RUN and ANALYZE MODIS L1A Data ----


clear variables;
scriptPlotting_blk;
% By Andrew J. Buggee
%% ---- Read Data ----

folderName = './MODIS_data/2017_5_14/L1B_calibrated_radiances/';
fileName = 'MOD02HKM.A2021237.1920.006.2021238072711.hdf';
geoFileName = 'MOD03.A2021237.1920.006.2021238012056.hdf';

[EV,geo,solar,sensor] = readMODIS_L1A_data(folderName,fileName,geoFileName);



%% plot frames


for ii = 1:size(EV.m250,3)
    figure;
    subplot(1,2,1)
    imagesc(EV.m250(:,:,ii));
    colorbar
    title(['Band: ',num2str(sensor.modisBands(ii)),' \mum'])
    pause(0.33)
    
    % plot the histogram of counts
    nbins = 100;
    subplot(1,2,2)
    histogram(EV.m250(:,:,ii),nbins)
    xlabel('Raw Counts')
    ylabel('Frequency')
    set(gca,'yscale','log')
end


for ii = 1:size(EV.m500,3)
    
    figure;
    subplot(1,2,1)
    imagesc(EV.m500(:,:,ii));
    colorbar
    title(['Band: ',num2str(sensor.modisBands(ii+2)),' \mum'])
    pause(0.33)
    
        % plot the histogram of counts
    nbins = 100;
    subplot(1,2,2)
    histogram(EV.m500(:,:,ii),nbins)
    xlabel('Raw Counts')
    ylabel('Frequency')
    set(gca,'yscale','log')
end