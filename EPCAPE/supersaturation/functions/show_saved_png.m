function fig = show_saved_png(png_path, width_px)
% SHOW_SAVED_PNG  Show a saved figure PNG so that it fills its window, with no margins.
%
%   fig = show_saved_png(png_path)
%   fig = show_saved_png(png_path, width_px)
%
%   The figures in the live script are drawn invisibly and saved (exporting a visible
%   figure has hung MATLAB R2024b), and the saved PNG is shown instead. Shown with a
%   plain image(imread(...)), the image sits in a default axes, whose position
%   [0.13 0.11 0.775 0.815] leaves margins of 11-13 % on every side, and "axis image"
%   then shrinks it further whenever the window's shape differs from the image's
%   (a 2.1 : 1 image in a ~1.3 : 1 window loses ~40 % of the height). Here the window
%   gets the image's aspect ratio and the axes fill it, so the image is as large as
%   the window allows.
%
%   Inputs
%     png_path  the PNG to show
%     width_px  window width in pixels (default 1250); the height follows from the
%               image's aspect ratio
%
%   One image object, so it is cheap to display in the Live Editor.

if nargin < 2 || isempty(width_px)
    width_px = 1250;
end

img = imread(png_path);
[h_img, w_img, ~] = size(img);                        % image size, pixels
height_px = round(width_px * h_img / w_img);          % window height with the image's aspect ratio

fig = new_figure('Position', [100 100 width_px height_px]);
ax = axes(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);   % axes fill the whole window
image(ax, img);
axis(ax, 'image', 'off');

end
