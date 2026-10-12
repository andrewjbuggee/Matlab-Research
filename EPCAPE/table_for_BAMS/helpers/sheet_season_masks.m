function masks = sheet_season_masks(t)
% SHEET_SEASON_MASKS  Which samples fall in each of the table's five columns.
%
%   masks = sheet_season_masks(t)
%
%   Input
%     t      datetime vector (UTC) of the samples
%   Output
%     masks  logical matrix, one row per sample and one column per season,
%            in the order of the table's columns E-I:
%              column 1  EPCAPE  2023-02-15 .. 2024-02-14   (whole campaign)
%              column 2  Spring  2023-04-01 .. 2023-06-30   (header "AMJn23")
%              column 3  Summer  2023-07-01 .. 2023-09-30   (header "JlAS23")
%              column 4  Fall    2023-10-01 .. 2023-12-31   (header "OND23")
%              column 5  Winter  2023-02-15 .. 2023-03-31  AND  2024-01-01 .. 2024-02-14
%                                (header "FM23+JF24")
%
%   These are the seasons written in the sheet's column headers. They are
%   applied in UTC. Every day is included whole: "2023-06-30" means up to,
%   but not including, 2023-07-01 00:00.
%
%   Some of Kavin's own calculations use calendar quarters instead (winter =
%   Jan-Mar); see calendar_quarter_masks.m.

t = t(:);   % make sure t is a column

% A helper that is true for times in [first_day 00:00, last_day + 1 day 00:00)
in_days = @(first_day, last_day) (t >= datetime(first_day)) & (t < datetime(last_day) + days(1));

masks = false(numel(t), 5);
masks(:, 1) = in_days('2023-02-15', '2024-02-14');                                   % EPCAPE
masks(:, 2) = in_days('2023-04-01', '2023-06-30');                                   % Spring
masks(:, 3) = in_days('2023-07-01', '2023-09-30');                                   % Summer
masks(:, 4) = in_days('2023-10-01', '2023-12-31');                                   % Fall
masks(:, 5) = in_days('2023-02-15', '2023-03-31') | in_days('2024-01-01', '2024-02-14');  % Winter

end
