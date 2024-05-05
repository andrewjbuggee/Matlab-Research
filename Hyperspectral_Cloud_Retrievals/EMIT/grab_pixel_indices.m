%% Grab the pixel indices for each row and column provided


% By Andrew John Buggee

function pixels2use = grab_pixel_indices(pixels2use, matrix_size)

if length(pixels2use.row)==1 && length(pixels2use.col)==1

    pixels2use.idx = sub2ind(matrix_size, pixels2use.row, pixels2use.col);    

else

    for rr = 1:length(pixels2use.row)
        for cc = 1:length(pixels2use.col)

            pixels2use.idx(rr*cc) = sub2ind(matrix_size, pixels2use.row(rr), pixels2use.col(cc));    

        end
    end
    


end
