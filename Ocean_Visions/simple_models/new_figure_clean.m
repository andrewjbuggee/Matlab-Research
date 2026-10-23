function fig = new_figure_clean(varargin)
% NEW_FIGURE_CLEAN  An INVISIBLE figure with MATLAB's factory graphics defaults.
%
%   fig = new_figure_clean('Position', [100 100 1200 700])
%
% Two things this does, both learned the hard way (see the EPCAPE BAMS Fig. 7
% and supersaturation work, Oct 2026; template: EPCAPE/supersaturation/
% functions/new_figure.m):
%
%   1. startup.m sets dark ROOT defaults (black background, white 18-pt bold
%      text, 4-pt lines). Instead of changing the root for the whole session,
%      every default the root carries is reset to its FACTORY value on this one
%      figure, so everything drawn in it afterwards starts from the factory
%      look and other figures are unaffected.
%   2. The figure is created 'Visible','off'. R2024b stalls (idle CPU, >10 min)
%      when exporting VISIBLE figures from -batch runs; invisible figures export
%      in seconds. Save with exportgraphics, then close(fig).

fig = figure('Visible', 'off', varargin{:});
defs = fieldnames(get(groot, 'default'));
for k = 1:numel(defs)
    name = defs{k};
    if startsWith(name, 'defaultFigure', 'IgnoreCase', true)
        continue   % figure-level properties are set explicitly below
    end
    set(fig, name, get(groot, regexprep(name, '^default', 'factory', 'ignorecase')));
end
fig.Color = 'w';
end
