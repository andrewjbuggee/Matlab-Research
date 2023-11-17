% grab modis reflectacne data for the bands of interest and the pixels
% being used


% By Andrew J. Buggee
%%

function modisRefl = grab_modis_reflectance(modis,inputs, pixels2use)

% define the MODIS bands to use
bands2run = inputs.bands2run;

% define the number of pixels used
numPixels = inputs.pixels.num_2calculate;

% define the zero array of MODIS reflectance
modisRefl = zeros(numPixels, length(bands2run));


for pp = 1:numPixels

    for bb = 1:length(bands2run)
        
        modisRefl(pp,bb) = modis.EV1km.reflectance(pixels2use.res1km.row(pp),pixels2use.res1km.col(pp),bands2run(bb));

    end

    
end




end