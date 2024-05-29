%% Grab the pixel indices for each row and column provided


% By Andrew John Buggee

function pixels2use = grab_pixel_indices(pixels2use, matrix_size)



if length(pixels2use.row)==1 && length(pixels2use.col)==1

    pixels2use.idx = sub2ind(matrix_size, pixels2use.row, pixels2use.col);    

else

    pixels2use.idx = zeros(1, length(pixels2use.row));

    for ii = 1:length(pixels2use.idx)

            pixels2use.idx(ii) = sub2ind(matrix_size, pixels2use.row(ii), pixels2use.col(ii));    

    end
    


end
