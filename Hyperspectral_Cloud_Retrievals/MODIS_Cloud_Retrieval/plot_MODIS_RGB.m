function [f] = plot_MODIS_RGB(modis, inputs)


R = modis.EV1km.radiance(:,:,2);      % 850 nm band is going to be our red hue
G = modis.EV1km.radiance(:,:,1);      % 650 nm band is going to be our green hue
B = modis.EV1km.radiance(:,:,3);      % 465 nm band is going to be our blue hue

numRows = size(R,1);
numCols = size(R,2);

% If there are any radiance values less than 0, set them to be 0
% R(R<0) = 0;
% G(G<0) = 0;
% B(B<0) = 0;

% Lets normalize each band and convert it to 8-bit image

R = uint8(R);
G = uint8(G);
B = uint8(B);

% create the iamge array
I = cat(3,R,G,B);

% for some reason we our aray is upside down. We want the equatorward
% direction to be south, and the polewrad direction to be north

I = rot90(I,2);

% lets make a plot!

f = figure; 
image(I); 
title(['True Color'])
set(f, 'Position', [0 0 numCols numRows])

end

