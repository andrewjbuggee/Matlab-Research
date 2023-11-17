%% ---- Get MODIS .INP Names based on pixel -----



% By Andrew J. Buggee

%%

function [inpNames] = getMODIS_INPnames_noClouds(solar,pixel_row,pixel_column,bands2run)

rowName = num2str(pixel_row);
colName = num2str(pixel_column);

szaName = num2str(solar.zenith(pixel_row,pixel_column)./100);
sazName = num2str(solar.azimuth(pixel_row,pixel_column)./100);

inpNames = cell(1,length(bands2run));

for ii = 1:length(bands2run)
    
    bandName = num2str(bands2run(ii));
    inpNames{ii} = ['pixel_',rowName,'r_',colName,...
        'c_sza_',szaName,'_saz_',sazName,'_band_',bandName,'.INP'];
    
end





end
