function f = plot_R16_and_R17_withUnertainty(modis,inputs)


liquidWater_mask = modis.cloud.phase == 2; % 2 is the value designated for liquid water

% create tau mask based on threshold
tauThreshold = inputs.pixels.tauThreshold;

tau_mask = modis.cloud.optThickness16 >= tauThreshold;




index_re = modis.cloud.effRadius16>=0 & modis.cloud.effRadius17>=0;        % find values greater than 0

% find where there is overlap

combined_mask = logical(liquidWater_mask .* tau_mask.*index_re);

re16 = modis.cloud.effRadius16(combined_mask);

re17 = modis.cloud.effRadius17(combined_mask);

yneg = re17.*(modis.cloud.effRad_uncert_17(combined_mask)./100);
ypos = yneg;
xneg = re16.*(modis.cloud.effRad_uncert_16(combined_mask)./100);
xpos = xneg;

% ----- THERE ARE TWO MANY PIXELS!! FIND A RANDOM SAMPLE -----

num2plot = 1e3;

if numel(re16)>num2plot
    
    
    rand_ind = randi([1,numel(re16)],num2plot,1);
    
    f = figure; errorbar(re16(rand_ind),re17(rand_ind),yneg(rand_ind),ypos(rand_ind),xneg(rand_ind),xpos(rand_ind),'o')
    xlabel('r_{e}^{16} (\mu m)')
    ylabel('r_{e}^{17} (\mu m)')
    title('MODIS Droplet Radius Estimates')
    grid on; grid minor
    set(f, 'Position', [0 0 1000 400])
    
else
    f = figure; errorbar(re16,re17,yneg,ypos,xneg,xpos,'o')
    xlabel('r_{e}^{16} (\mu m)')
    ylabel('r_{e}^{17} (\mu m)')
    title('MODIS Droplet Radius Estimates')
    grid on; grid minor
    set(f, 'Position', [0 0 1000 400])
    
end



end

