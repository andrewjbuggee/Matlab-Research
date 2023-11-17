%% ----- Extract MODIS calculated values as truth -----



% By Andrew J. Buggee

%%

function [truth_estimate_table] = gatherTruthEstimateVals(modis, minVals, inputs,pixels2use)

% extract inputs

bands2run = inputs.bands2run;
bands2plot = inputs.bands2plot;
num_pixels = inputs.pixels.num_2calculate;
bands2search = inputs.bands2search;

% 1km resolution pixels
pixel_row_1 = pixels2use.res1km.row;
pixel_col_1 = pixels2use.res1km.col;

% 500 meter resolution pixels
%pixel_row_500 = pixels2use.res500m.row;
%pixel_col_500 = pixels2use.res500m.col;

% save calculations
saveCalcs_filename = inputs.saveCalculations_fileName;


% create reflectance function table
truth_estimate_table = table;



% extract the modis computed values for reflectance to compare with my
% internal computation of reflectance

% SWITCH TO 1KM DATA
% but for now, lets propose a simple algebraic fix so taht we select
% the proper pixel in our 500 meter data set
% lets take an average of the four pixels that map to the 500 meter
% data
pixelIndex_500m = 1:4:((num_pixels-1)*4 +1);

for pp = 1:num_pixels
   

    
    
    
    
    % extract the values MODIS calculates using its own algorithm
    truth_estimate_table.modisR17(pp) = modis.cloud.effRadius17(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisR17_uncert(pp) = modis.cloud.effRad_uncert_17(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisT17(pp) = modis.cloud.optThickness17(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisT17_uncert(pp) = modis.cloud.optThickness_uncert_17(pixel_row_1(pp),pixel_col_1(pp));
    
    truth_estimate_table.modisR16(pp) = modis.cloud.effRadius16(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisR16_uncert(pp) = modis.cloud.effRad_uncert_16(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisT16(pp) = modis.cloud.optThickness16(pixel_row_1(pp),pixel_col_1(pp));
    truth_estimate_table.modisT16_uncert(pp) = modis.cloud.optThickness_uncert_16(pixel_row_1(pp),pixel_col_1(pp));
    
    % gather the estimated values in the table
    
    truth_estimate_table.estR17(pp) = minVals.minR(1,pp);  % My estimates for bands 1 and 7
    truth_estimate_table.estT17(pp) = minVals.minT(1,pp);
    
%     truth_estimate_table.estR27(pp) = minVals.minR(2,pp); % My estiamtes for bands 2 and 7 (should have a different tau)
%     truth_estimate_table.estT27(pp) = minVals.minT(2,pp);
    
%     truth_estimate_table.estR16(pp) = minVals.minR(2,pp); % My estimates for bands 1 and 6 ( should have a different estiamte for re)
%     truth_estimate_table.estT16(pp) = minVals.minT(2,pp);
    
    % compute the absolute difference
    truth_estimate_table.squareDiffR17(pp) = (minVals.minR(1,pp) - truth_estimate_table.modisR17(pp)).^2;
    truth_estimate_table.squareDiffT17(pp) = (minVals.minT(1,pp) - truth_estimate_table.modisT17(pp)).^2;
    
    % compute the absolute difference
%     truth_estimate_table.squareDiffR27(pp) = (minVals.minR(2,pp) - truth_estimate_table.modisR17(pp)).^2;
%     truth_estimate_table.squareDiffT27(pp) = (minVals.minT(2,pp) - truth_estimate_table.modisT17(pp)).^2;
    
    % compute the absolute difference
%     truth_estimate_table.squareDiffR16(pp) = (minVals.minR(2,pp) - truth_estimate_table.modisR16(pp)).^2;
%     truth_estimate_table.squareDiffT16(pp) = (minVals.minT(2,pp) - truth_estimate_table.modisT16(pp)).^2;
    
    
end





save(saveCalcs_filename,"truth_estimate_table",'-append'); % save inputSettings to the same folder as the input and output file






end
