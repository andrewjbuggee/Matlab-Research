%% Create a figure that uses MATLAB's factory graphics defaults
%
% DESCRIPTION:
%   startup.m sets dark root graphics defaults (black background, white 18-pt
%   bold text, line width 4). Those defaults make the multi-panel adiabaticity
%   figures unreadable: axis labels and titles collide, and light series
%   colors lose contrast. Rather than changing the root defaults for the whole
%   session, each figure gets the factory values as its OWN defaults, so every
%   object drawn in it afterward starts from the factory values and all other
%   figures are unaffected.
%
%   Same recipe as EPCAPE/supersaturation/functions/new_figure.m. Create the
%   figure with this function BEFORE drawing in it.
%
% INPUTS:
%   varargin - anything figure() accepts, e.g. new_factory_figure('Visible', 'off')
%
% OUTPUT:
%   fig - the figure handle
%
% By Andrew John Buggee
%%

function fig = new_factory_figure(varargin)

fig = figure(varargin{:});

% The defaults that startup.m (or anything else) has set at the root level
root_defaults = fieldnames(get(groot, 'default'));

for kk = 1:numel(root_defaults)

    name = root_defaults{kk};

    % Figure-level properties belong to the figure itself and are set below
    if startsWith(name, 'defaultFigure', 'IgnoreCase', true)
        continue
    end

    set(fig, name, get(groot, regexprep(name, '^default', 'factory', 'ignorecase')));

end

fig.Color = 'w';

end
