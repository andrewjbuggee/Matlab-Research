%% ---- Get MODIS .INP Names based on pixel -----



% By Andrew J. Buggee

%%

function [inpNames] = getMODIS_INPnames_withClouds(solar,inputs,pixels2use)

% extract inputs

num_pixels = inputs.pixels.num_2calculate;
re = inputs.re;
tau_c = inputs.tau_c;
bands2run = inputs.bands2run;
pixel_row = pixels2use.res1km.row;
pixel_col = pixels2use.res1km.col;






inpNames = cell(num_pixels,length(re),length(tau_c),length(bands2run));



for pp = 1:num_pixels
    
    
    rowName = num2str(pixel_row(pp));
    colName = num2str(pixel_col(pp));
    
    szaName = num2str(solar.zenith(pixel_row(pp),pixel_col(pp)));
    sazName = num2str(solar.azimuth(pixel_row(pp),pixel_col(pp)));
    
    for ii = 1:length(bands2run)
        
        for jj = 1:length(re)
            
            for kk = 1:length(tau_c)
                
                bandName = num2str(bands2run(ii));
                inpNames{pp,jj,kk,ii} = ['pixel_',rowName,'r_',colName,...
                    'c_sza_',szaName,'_saz_',sazName,'_band_',bandName,...
                    '_r_',num2str(re(jj)),'_T_',num2str(tau_c(kk)),'.INP'];
                
            end
            
        end
        
    end
    
end



end
