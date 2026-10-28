%% Which minutes belong to a chosen set of cloud cases (Dedrick's 24 events, Soledad in-cloud, low-cloud periods, ...)
%
% DESCRIPTION:
%   Returns a logical mask on the 1-min grid of build_epcape_1min_inputs.m
%   for one of these case sets:
%
%   "dedrick24"        the in-cloud windows of the 24 Mt. Soledad cloud
%                      events of Dedrick et al. (2026, supplement Table S5),
%                      via filter_to_dedrick2026_events.m. With
%                      cloudIntervalsOnly = true, only the 5-min intervals of
%                      each window that meet the cloud condition are kept.
%   "soledadInCloud"   every minute in a 5-min interval meeting Dedrick's
%                      cloud condition at Mt. Soledad (visibility < 1 km,
%                      FM-120 LWC > 0.01 g m-3), event or not: the 24
%                      events with ALL of their GCVI/SEMS/duration criteria
%                      relaxed
%   "lowCloudPeriods"  Dedrick's hourly low-cloud periods at the Pier
%                      (single layer, base and top < 3 km, |base - LCL| <
%                      150 m, i.e. coupled; cloud NOT enveloping Mt.
%                      Soledad; supplement Text S7, Table S6)
%   "all"              every minute (only the analysis filters select)
%
% INPUTS:
%   D                  - 1-min timetable (build_epcape_1min_inputs.m)
%   caseSet            - one of the strings above
%   cloud5             - 5-min Mt. Soledad classification (for cloudIntervalsOnly)
%   H                  - hourly timetable from load_hourly_ccn_and_supersaturation.m
%                        (for its low_cloud_period flag)
%   cloudIntervalsOnly - (optional, default false) see "dedrick24"
%
% OUTPUTS:
%   in_case  - logical column, one per minute
%   event_id - Dedrick event number of each minute (NaN outside; only for
%              "dedrick24", NaN otherwise)
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Text S7, Tables S4-S6.
%
% By Andrew John Buggee
%%

function [in_case, event_id] = case_set_mask_1min(D, caseSet, cloud5, H, cloudIntervalsOnly)


if nargin < 5 || isempty(cloudIntervalsOnly)
    cloudIntervalsOnly = false;
end
n = height(D);
event_id = nan(n, 1);

switch lower(string(caseSet))

    case "dedrick24"
        % each minute [t, t + 1 min) is kept if >= half of it lies in a window
        args = {'Period', "inCloud", 'SampleLength', minutes(1), 'TimeLabel', "start", 'MinOverlap', 0.5};
        if cloudIntervalsOnly
            args = [args, {'SoledadCloud5min', cloud5}];
        end
        [in_case, event_id] = filter_to_dedrick2026_events(D.time, args{:});

    case "soledadincloud"
        in_case = D.soledad_cloud == 1;

    case "lowcloudperiods"
        % the hour [h, h + 1) of each minute
        h = dateshift(D.time, 'start', 'hour');
        [tf, loc] = ismember(h, H.time);
        in_case = false(n, 1);
        in_case(tf) = H.low_cloud_period(loc(tf)) == 1;

    case "all"
        in_case = true(n, 1);

    otherwise
        error('case_set_mask_1min: unknown case set "%s"', caseSet);
end


end
