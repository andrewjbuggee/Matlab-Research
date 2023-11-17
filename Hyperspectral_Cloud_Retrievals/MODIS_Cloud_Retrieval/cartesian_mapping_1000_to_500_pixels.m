% ---- MAPPING between two pixel grids -----



% By Andrew J. Buggee

%%

function [newRows,newCols] = cartesian_mapping_1000_to_500_pixels(pixels)


row1 = pixels.res1km.size(1);
col1 = pixels.res1km.size(2);

row2 = pixels.res500m.size(1);
col2 = pixels.res500m.size(2);

rowFactor = row2/row1;
colFactor = col2/col1;


numPixels = length(pixels.res1km.col);

newIndex = [];

for ii = 1:numPixels
    

    % if our larger pixel is in the ith row, we add i-1 to find the row for
    % the smaller pixels
    
    addRow = pixels.res1km.row(ii) + [pixels.res1km.row(ii)-1, pixels.res1km.row(ii)];
    addCol = pixels.res1km.col(ii) + [pixels.res1km.col(ii)-1, pixels.res1km.col(ii)];
    

    newIndex = [newIndex;[addRow(1),addCol(1);addRow(2),addCol(1); addRow(1),addCol(2);addRow(2),addCol(2)]];




end


newRows = newIndex(:,1); 

newCols = newIndex(:,2);




end
