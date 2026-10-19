function fig = new_figure(varargin)
% NEW_FIGURE  A figure that uses MATLAB's FACTORY graphics defaults, whatever startup.m sets.
%
%   fig = new_figure(...)        same inputs as figure(), e.g. new_figure('Visible', 'off', 'Position', [...])
%
%   Andrew's startup.m sets dark root defaults (black background, white text,
%   18-pt bold fonts, line width 4). The plots here use the paper's light and
%   dark greys, which need a white background. So, instead of changing the root
%   defaults for the whole session, each figure gets the factory values as its
%   OWN defaults: every object drawn in it afterwards (axes, lines, text,
%   legends) starts from the factory values, and other figures are unaffected.
%   (Recipe from the EPCAPE BAMS Fig. 7 work, Oct 2026.)
%
%   Create every figure with this function BEFORE drawing in it.

fig = figure(varargin{:});
defs = fieldnames(get(groot, 'default'));   % the defaults startup.m (or anything else) has set
for k = 1:numel(defs)
    name = defs{k};
    if startsWith(name, 'defaultFigure', 'IgnoreCase', true)
        continue   % figure-level properties of the figure itself are set below
    end
    set(fig, name, get(groot, regexprep(name, '^default', 'factory', 'ignorecase')));
end
fig.Color = 'w';

end
