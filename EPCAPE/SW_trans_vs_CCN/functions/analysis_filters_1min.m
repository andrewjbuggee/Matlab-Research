%% Evaluate the cloud-scene filters minute by minute (daytime, LWP, no rain, single layer, liquid, persistent)
%
% DESCRIPTION:
%   Takes the 1-min grid of build_epcape_1min_inputs.m and the thresholds
%   set at the top of the analysis scripts, and returns one logical mask per
%   filter (true = the minute passes). Every filter can be switched off.
%   The masks are kept separate so the scripts can count how much data each
%   one removes, alone and in sequence.
%
%   THE FILTERS (field names of the output in brackets)
%
%   [daytime] the sun is higher than 90 - szaMax_deg degrees:
%       cos(SZA) > cos(szaMax_deg). The SW transmittance needs sunlight, and
%       at low sun the ratio is noisy (sw_transmittance.m).
%
%   [lwp] lwpMin_gm2 < LWP < lwpMax_gm2 (MWRRET v2, cloudy samples, QC).
%       The lower bound is about the retrieval's uncertainty (20 g m-2 is
%       the threshold the ARM NDROP VAP uses, "liquid water path above
%       measurement uncertainty"; Dedrick et al. 2026, supplement Text S5).
%       The upper bound removes clouds likely to precipitate (NDROP uses
%       300 g m-2; 250 is a stricter choice). A minute with no LWP sample
%       fails.
%
%   [noPrecip] no rain at the chosen disdrometers within +/- rainWindow_min
%       of the minute. The window also removes the time a wet radome
%       inflates the microwave LWP (23.8 and 31.4 GHz emission from liquid
%       on the radome; the +/- 30 min screen of quantities.near_rain).
%       Minutes without disdrometer data are NOT removed (no evidence of
%       rain); the scripts report how often that happens.
%
%   [singleLayer] no evidence of a second cloud layer: KAZR-ARSCL reports
%       at most one hydrometeor layer in every 4-s sample of the minute, and
%       the ceilometer reports no second or third base (detection_status
%       2-3). This combines Dedrick's Method 2 (ARSCL: only one layer in the
%       column) with the ceilometer, because ARSCL's KAZR-based layers miss
%       ~1/3 of the lidar's low cloud bases at EPCAPE (EPCAPE CloudCube work,
%       Oct 2026), so "exactly one ARSCL layer" alone would also drop real
%       single-layer stratocumulus. If requireArscl is true a minute without
%       ARSCL data fails. NOTE: Dedrick used Method 1, the NDROP VAP's
%       cloud_base_type = 1, which is not in the local data.
%
%   [liquid] a PROXY for liquid phase: the highest ARSCL layer top <=
%       maxCloudTop_m and the ceilometer cloud base <= maxCloudBase_m
%       (Dedrick's low-cloud criterion, "cloud top and base heights < 3 km",
%       supplement Text S7). No cloud-top temperature is available locally
%       (INTERPSONDE is not downloaded; only 6 radiosondes are), so phase
%       cannot be tested directly. A 3-km top is below the 0 deg C level of
%       the coastal southern California boundary layer for most of the year
%       (my inference, not checked against soundings), so these clouds are
%       warm; check winter cases with soundings before trusting this proxy.
%       A minute with no ARSCL layer passes the top test (no evidence of a
%       high top) but needs a ceilometer base.
%
%   [persistence] a cloud deck is present at each site in persistenceSites
%       for >= persistMinFraction of the window minute +/- persistHalfWindow_min:
%         Pier:    a minute is cloudy if >= 50 % of its ceilometer samples
%                  see a cloud base (pier_cloud_frac >= 0.5)
%         Soledad: a minute is cloudy if it lies in a 5-min interval meeting
%                  Dedrick's cloud condition (visibility < 1 km throughout,
%                  FM-120 LWC > 0.01 g m-3). This means Mt. Soledad is INSIDE
%                  the cloud (base below ~250 m MSL): with the local data, a
%                  deck whose base is above Mt. Soledad cannot be seen there.
%       The fraction is over the window minutes that have data, and at
%       least persistMinCoverage of the window must have data.
%
% INPUTS:
%   D - timetable from build_epcape_1min_inputs.m
%   F - struct of settings (all fields below are needed):
%       applyDaytime, szaMax_deg
%       applyLwp, lwpMin_gm2, lwpMax_gm2
%       applyNoPrecip, rainWindow_min, rainSites ("Pier", "Soledad")
%       applySingleLayer, requireArscl
%       applyLiquid, maxCloudTop_m, maxCloudBase_m
%       applyPersistence, persistenceSites ("Pier", "Soledad"),
%           persistHalfWindow_min, persistMinFraction, persistMinCoverage
%
% OUTPUTS:
%   M - struct of logical column vectors (one element per minute): daytime,
%       lwp, noPrecip, singleLayer, liquid, persistence, and all (the AND of
%       the filters that are switched on). A switched-off filter is all true.
%       M.applied lists the switched-on filters, M.description says in words
%       what each one tested.
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1; supplement Text S5, S7, Tables S4, S6.
%   Zhang, D., and coauthors, 2023: Evaluation of four ground-based
%     retrievals of cloud droplet number concentration in marine
%     stratocumulus with aircraft in situ measurements. Atmos. Meas. Tech.,
%     16, 5827-5846, doi:10.5194/amt-16-5827-2023 (the NDROP VAP screens).
%
% By Andrew John Buggee
%%

function M = analysis_filters_1min(D, F)


n = height(D);
all_true = true(n, 1);
M = struct();
M.description = struct();



% -------------------------------------------------------------------
% ------------------------------ Daytime ----------------------------
% -------------------------------------------------------------------

M.daytime = all_true;
if F.applyDaytime
    M.daytime = D.cos_sza > cosd(F.szaMax_deg);
end
M.description.daytime = sprintf('solar zenith angle < %g deg', F.szaMax_deg);



% -------------------------------------------------------------------
% ------------------------ Liquid water path ------------------------
% -------------------------------------------------------------------

M.lwp = all_true;
if F.applyLwp
    M.lwp = D.lwp_gm2 > F.lwpMin_gm2 & D.lwp_gm2 < F.lwpMax_gm2;   % NaN (no sample) fails
end
M.description.lwp = sprintf('%g < LWP < %g g m-2 (MWRRET v2, cloudy, QC)', F.lwpMin_gm2, F.lwpMax_gm2);



% -------------------------------------------------------------------
% ---------------------- No precipitation nearby --------------------
% -------------------------------------------------------------------

M.noPrecip = all_true;
if F.applyNoPrecip
    rain_any = false(n, 1);
    for site = string(F.rainSites)
        switch lower(site)
            case "pier",    r = D.rain_pier;
            case "soledad", r = D.rain_soledad;
            otherwise, error('analysis_filters_1min: unknown rain site "%s"', site);
        end
        rain_any = rain_any | (r == 1);              % NaN (no data) is not rain
    end
    % rain anywhere in [t - w, t + w] removes minute t
    near = movmax(double(rain_any), [F.rainWindow_min F.rainWindow_min]) > 0;
    M.noPrecip = ~near;
end
M.description.noPrecip = sprintf('no rain at %s within +/- %g min', strjoin(string(F.rainSites), ' & '), F.rainWindow_min);



% -------------------------------------------------------------------
% --------------------------- Single layer --------------------------
% -------------------------------------------------------------------

M.singleLayer = all_true;
if F.applySingleLayer
    arscl_ok = D.arscl_max_layers <= 1;                     % NaN (no ARSCL) -> false
    if ~F.requireArscl
        arscl_ok = arscl_ok | isnan(D.arscl_max_layers);    % no ARSCL data: not evidence of a 2nd layer
    end
    ceil_ok = ~(D.pier_multi_base == 1);                    % NaN (no ceilometer) -> no evidence
    M.singleLayer = arscl_ok & ceil_ok;
end
M.description.singleLayer = 'ARSCL <= 1 layer and no 2nd/3rd ceilometer base';



% -------------------------------------------------------------------
% ----------------------- Liquid (height proxy) ---------------------
% -------------------------------------------------------------------

M.liquid = all_true;
if F.applyLiquid
    top_ok  = D.arscl_top_max_m <= F.maxCloudTop_m | isnan(D.arscl_top_max_m);
    base_ok = D.pier_cbh_m <= F.maxCloudBase_m;             % NaN (no base) fails
    M.liquid = top_ok & base_ok;
end
M.description.liquid = sprintf('proxy: ARSCL top <= %g m and ceilometer base <= %g m', F.maxCloudTop_m, F.maxCloudBase_m);



% -------------------------------------------------------------------
% ---------------- Persistence of the deck at each site -------------
% -------------------------------------------------------------------

M.persistence = all_true;
if F.applyPersistence
    w = F.persistHalfWindow_min;
    n_window = 2 * w + 1;
    for site = string(F.persistenceSites)
        switch lower(site)
            case "pier"
                c = double(D.pier_cloud_frac >= 0.5);
                c(isnan(D.pier_cloud_frac)) = NaN;
            case "soledad"
                c = D.soledad_cloud;                         % 1 / 0 / NaN
            otherwise
                error('analysis_filters_1min: unknown persistence site "%s"', site);
        end
        has = ~isnan(c);
        n_has   = movsum(double(has), [w w]);                % window minutes with data
        n_cloud = movsum(double(c == 1), [w w]);             % window minutes with cloud
        frac = n_cloud ./ n_has;
        coverage = n_has / n_window;
        M.persistence = M.persistence & frac >= F.persistMinFraction & coverage >= F.persistMinCoverage;
    end
end
M.description.persistence = sprintf('cloud at %s for >= %g %% of +/- %g min', ...
    strjoin(string(F.persistenceSites), ' & '), 100 * F.persistMinFraction, F.persistHalfWindow_min);



% -------------------------------------------------------------------
% ------------------------- All filters together --------------------
% -------------------------------------------------------------------

names = ["daytime" "lwp" "noPrecip" "singleLayer" "liquid" "persistence"];
on    = [F.applyDaytime F.applyLwp F.applyNoPrecip F.applySingleLayer F.applyLiquid F.applyPersistence];
M.applied = names(on);
M.all = all_true;
for name = M.applied
    M.all = M.all & M.(name);
end


end
