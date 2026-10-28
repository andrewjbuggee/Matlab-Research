%% Hourly SW transmittance, LWP and sun angle from the minutes that pass every filter
%
% DESCRIPTION:
%   The supersaturation retrievals are hourly (Dedrick et al. 2026 averaged
%   every retrieval to 1 h), and the Pier CCN spectra are assembled hourly,
%   so the scatter plot is hourly. This function turns the 1-min grid into
%   hours [h, h + 1) (UTC, labelled by the start):
%     * T_sw: sw_transmittance.m (SKYRAD SW / RADFLUX clear-sky SW, ratio of
%       the hour's means) over ONLY the minutes that pass, so cloud-free or
%       rainy minutes inside an otherwise good hour do not enter the ratio
%     * lwp_gm2: mean MWRRET v2 LWP of the same minutes
%     * sza_deg: solar zenith angle of the mean cos(SZA) of those minutes
%     * n_pass: how many minutes passed
%   An hour is kept (finite T_sw) only if at least min_pass_minutes of its
%   60 minutes pass every filter AND have valid SW and clear-sky values.
%
% INPUTS:
%   D                - 1-min timetable (build_epcape_1min_inputs.m)
%   pass             - logical, one per minute: in the case set AND passing
%                      every filter
%   min_pass_minutes - fewest passing minutes for an hour to count (e.g. 30)
%   szaMax_deg       - the daytime cut handed to sw_transmittance.m
%
% OUTPUTS:
%   Hr - hourly timetable: T_sw, lwp_gm2, sza_deg, n_pass
%
%   Kept for dedrick2026_data_loss_analysis.m; sw_transmittance_vs_ccn.m
%   calls average_passing_minutes.m with any window length.
%
% DEPENDENCIES: average_passing_minutes.m, EPCAPE/SW_transmittance
%
% By Andrew John Buggee
%%

function Hr = hourly_from_passing_minutes(D, pass, min_pass_minutes, szaMax_deg)


% one-hour windows of average_passing_minutes.m (the clock hours [h, h + 1))
Hr = average_passing_minutes(D, pass, 60, min_pass_minutes, szaMax_deg);

end
