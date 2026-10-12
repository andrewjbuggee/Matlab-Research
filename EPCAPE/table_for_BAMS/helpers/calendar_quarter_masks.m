function masks = calendar_quarter_masks(t)
% CALENDAR_QUARTER_MASKS  Kavin's seasons: calendar quarters of any year.
%
%   masks = calendar_quarter_masks(t)
%
%   Several of Kavin's calculations (the SW, LW, FM-120 hours and the
%   15-min-file rows) do not use the sheet's seasons. They use the month
%   of each time stamp:
%     column 1  EPCAPE  every sample (no date limits at all)
%     column 2  Spring  April, May, June
%     column 3  Summer  July, August, September
%     column 4  Fall    October, November, December
%     column 5  Winter  January, February, March   (of BOTH years)
%
%   The only difference from sheet_season_masks.m is in winter and in the
%   campaign column: here winter is Jan-Mar (both years pooled, no date
%   limits), and the campaign column has no start or end date.
%
%   Input   t      datetime vector (the month is taken as given: pass UTC or
%                  local times, depending on what you want)
%   Output  masks  logical matrix (samples x 5), columns as listed above

m = month(t(:));   % calendar month (1-12) of every sample

masks = false(numel(m), 5);
masks(:, 1) = true;                       % EPCAPE: every sample
masks(:, 2) = ismember(m, [4 5 6]);       % Spring
masks(:, 3) = ismember(m, [7 8 9]);       % Summer
masks(:, 4) = ismember(m, [10 11 12]);    % Fall
masks(:, 5) = ismember(m, [1 2 3]);       % Winter

end
