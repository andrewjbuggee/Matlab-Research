%% Draw a filled, optionally stippled or striped, rectangle (patterns sized in inches on the page)
%
% DESCRIPTION:
%   MATLAB patches have no hatch fill, so the pattern is drawn as one line
%   object on top of a (semi-transparent) patch:
%     "none"     plain fill
%     "stipple"  dots on a staggered grid
%     "stripes"  45-degree stripes (45 deg on the page, whatever the axis
%                scales), clipped to the rectangle
%   Spacing is given in inches, so the pattern looks the same in every bar
%   regardless of the data units. This needs the axes' inches per data
%   unit, so the axes must have a fixed size (set its Position in inches
%   and its limits before calling this).
%
% INPUTS:
%   ax          - axes
%   x0, x1      - left and right edges (data units)
%   y0, y1      - bottom and top edges (data units)
%   color       - fill color [r g b]
%   alpha       - fill opacity (0-1)
%   pattern     - "none" | "stipple" | "stripes"
%   in_per_x    - inches per x data unit
%   in_per_y    - inches per y data unit
%   spacing_in  - pattern spacing (inches)
%   pat_color   - color of the dots / stripes
%
% OUTPUTS:
%   p - the patch (e.g. to add an outline)
%
% By Andrew John Buggee
%%

function p = draw_pattern_rect(ax, x0, x1, y0, y1, color, alpha, pattern, in_per_x, in_per_y, spacing_in, pat_color)

p = patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], color, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% rectangle in page inches (only the scale matters)
X0 = min(x0, x1) * in_per_x;  X1 = max(x0, x1) * in_per_x;
Y0 = min(y0, y1) * in_per_y;  Y1 = max(y0, y1) * in_per_y;
xs = [];  ys = [];

switch pattern
    case "none"
        return
    case "stipple"
        rows = (Y0 + spacing_in / 2):spacing_in:Y1;
        for r = 1:numel(rows)
            off = mod(r, 2) * spacing_in / 2;                     % stagger every other row
            cols = (X0 + spacing_in / 4 + off):spacing_in:X1;
            xs = [xs, cols]; %#ok<AGROW>
            ys = [ys, repmat(rows(r), 1, numel(cols))]; %#ok<AGROW>
        end
        if isempty(xs), return, end
        plot(ax, xs / in_per_x, ys / in_per_y, '.', 'Color', pat_color, 'MarkerSize', 2.5);
    case "stripes"
        for c = (X0 - (Y1 - Y0)):spacing_in:X1                  % lines X - (Y - Y0) = c
            xa = max(X0, c);  xb = min(X1, c + (Y1 - Y0));
            if xb <= xa, continue, end
            xs = [xs, xa, xb, NaN]; %#ok<AGROW>
            ys = [ys, Y0 + (xa - c), Y0 + (xb - c), NaN]; %#ok<AGROW>
        end
        if isempty(xs), return, end
        plot(ax, xs / in_per_x, ys / in_per_y, '-', 'Color', pat_color, 'LineWidth', 0.4);
    otherwise
        error('draw_pattern_rect: pattern must be "none", "stipple" or "stripes"');
end

end
