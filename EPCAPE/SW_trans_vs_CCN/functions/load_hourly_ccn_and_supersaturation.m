%% Hourly Pier CCN spectra and cloud supersaturations (Dedrick's posted values and my own s_Nd,FM)
%
% DESCRIPTION:
%   Builds one hourly timetable over the campaign with everything the
%   SW transmittance vs. CCN analysis needs about aerosol activation.
%
%   1. PIER CCN SPECTRA (cm-3) at the six CCN-counter supersaturations 0.1,
%      0.2, 0.4, 0.6, 0.8 and 1.0 %. Source: the DMT CCN-200 at the Scripps
%      Pier (ARM aosccn2colaavg), corrected by Dedrick et al. (2026) by a
%      factor 1.57 for an inlet leak that diluted the sample (supplement
%      Text S4), posted as epcccn2colaavgM1_DedrickCorrected.txt (local name
%      bb2743661p_35_1_CCN.txt), 30-min, online from 23 May 2023. This is
%      the CCN spectrum Dedrick et al. used for s_CCN and s_Nd. The counter
%      steps through the setpoints (~10 min each), so a 30-min row holds
%      1-3 of them; hourly_ccn_spectrum.m (supersaturation folder) averages
%      each setpoint over [h, h + 1), fills gaps of <= CcnGapFill_h hours by
%      linear interpolation in time, and keeps only complete spectra.
%      Uncertainty ~22 % on average, up to 60 % at 0.1-0.2 % (Text S3).
%
%   2. DEDRICK'S POSTED SUPERSATURATIONS (%), hourly, from
%      EPCAPE_Supersaturation.txt (bb2743661p_34_1_supersaturation.txt):
%      s_kappa, s_CCN, s_Nd and s_ARG00, each for "SoledadCloudEvent" and
%      "LowCloudPeriod" hours; NaN at every other hour. The skappa_* and
%      sCCN_* columns of the posted file are SWAPPED (evidence in
%      EPCAPE/supersaturation: s_CCN cannot exist before the CCN counter came
%      online on 23 May 2023, yet the posted sCCN_* columns start 12 Mar
%      2023; and the Fig. 10 panel statistics match the swapped labels).
%      With FixSwappedColumns = true (default) the labels are corrected here.
%
%   3. MY OWN s_Nd,FM ("AJB"), from the functions in EPCAPE/supersaturation
%      (reconstruction of Dedrick's cloud-scene method, r = 0.95 against his
%      posted s_Nd,FM; NOT yet independently verified): the hourly FM-120
%      droplet number at Mt. Soledad (mean of the 5-min values with LWC >
%      FmCloudLwc_gm3, i.e. in cloud) is located on the hour's Pier CCN
%      spectrum, extended linearly to (0 %, 0 cm-3) and to 1.2 %:
%      N_d = CCN(s_eff) (Twomey's activation argument; Rogers & Yau 1989,
%      pp. 95-96 and 107-110). It exists for every hour that Mt. Soledad is
%      in cloud and the Pier spectrum is complete, inside or outside the 24
%      events, which is what makes it useful when Dedrick's criteria are
%      relaxed. It is NOT available for low-cloud periods (cloud above Mt.
%      Soledad), which need the ARM NDROP droplet number (not downloaded).
%
% INPUTS:
%   products_folder - folder with Dedrick's posted products (bb2743661p_*)
%   fm              - FM-120 struct from read_fm120_daily_files.m
%   Optional name-value pairs:
%   'CcnGapFill_h'       2     (as in the supersaturation reconstruction)
%   'FmCloudLwc_gm3'     0.01  (FM-120 in-cloud threshold; Dedrick Text S7)
%   'FixSwappedColumns'  true
%   'TimeRange'          [datetime(2023,2,15) datetime(2024,2,15)]
%
% OUTPUTS:
%   H - hourly timetable (row time = hour start, UTC) with
%       ccn_cm3              (n x 6) complete Pier CCN spectrum, NaN rows otherwise
%       <method>_<period>    posted s (%), method = skappa, sCCN, sNd, sARG00;
%                            period = SoledadCloudEvent, LowCloudPeriod
%       low_cloud_period     Dedrick's hourly low-cloud flag (1/0)
%       nd_fm_cm3, n_fm_cloud   hourly in-cloud FM-120 N_d and the number of
%                            5-min intervals behind it
%       sNdFM_AJB            my s_Nd,FM (%)
%   s_set_pct - the six CCN setpoints (%), the columns of ccn_cm3
%
% DEPENDENCIES: EPCAPE/supersaturation/functions (read_posted_product,
%   hourly_ccn_spectrum, retrieve_supersaturation_pct).
%
% REFERENCES:
%   Dedrick, J. L., and coauthors, 2026: J. Atmos. Sci., 83, 1061-1082,
%     doi:10.1175/JAS-D-26-0004.1 (Sec. 2d, pp. 1065-1066) and supplement
%     Texts S3, S4, S7, S10.
%   Berta, V., S. Han, A. Osawa, J. Dedrick, N. Maneenoi, and L. M. Russell,
%     2025: Aerosol quantification products at Mt. Soledad and Scripps Pier
%     during EPCAPE 2023-24. UC San Diego Library Digital Collections,
%     doi:10.6075/J0WM1DSC (the posted CCN, supersaturation and low-cloud
%     files).
%   Rogers, R. R., and M. K. Yau, 1989: A Short Course in Cloud Physics,
%     3rd ed., Pergamon.
%
% By Andrew John Buggee
%%

function [H, s_set_pct] = load_hourly_ccn_and_supersaturation(products_folder, fm, varargin)


p = inputParser;
addParameter(p, 'CcnGapFill_h', 2);
addParameter(p, 'FmCloudLwc_gm3', 0.01);
addParameter(p, 'FixSwappedColumns', true);
addParameter(p, 'TimeRange', [datetime(2023, 2, 15) datetime(2024, 2, 15)]);
parse(p, varargin{:});
opt = p.Results;

s_set_pct = [0.1 0.2 0.4 0.6 0.8 1.0];   % Pier CCN-200 setpoints (%), the columns of the CCN file

% the hourly grid
t_hour = (opt.TimeRange(1):hours(1):opt.TimeRange(2) - hours(1))';
H = timetable(t_hour);
H.Properties.DimensionNames{1} = 'time';
n_h = numel(t_hour);
row_of_hour = @(t) round(hours(t - opt.TimeRange(1))) + 1;



% -------------------------------------------------------------------
% ---------------------- Pier CCN spectra (hourly) ------------------
% -------------------------------------------------------------------

CCN_30 = read_posted_product(fullfile(products_folder, 'bb2743661p_35_1_CCN.txt'));   % 30-min, cm-3
[t_ccn, ccn_h] = hourly_ccn_spectrum(CCN_30.time, CCN_30{:, :}, opt.CcnGapFill_h);
H.ccn_cm3 = nan(n_h, 6);
k = row_of_hour(t_ccn);
ok = k >= 1 & k <= n_h;
H.ccn_cm3(k(ok), :) = ccn_h(ok, :);



% -------------------------------------------------------------------
% ---------------- Dedrick's posted supersaturations ----------------
% -------------------------------------------------------------------

SS = read_posted_product(fullfile(products_folder, 'bb2743661p_34_1_supersaturation.txt'));   % hourly, %
if opt.FixSwappedColumns
    for period = ["SoledadCloudEvent" "LowCloudPeriod"]
        SS = renamevars(SS, "skappa_" + period, "TMP_" + period);
        SS = renamevars(SS, "sCCN_" + period, "skappa_" + period);
        SS = renamevars(SS, "TMP_" + period, "sCCN_" + period);
    end
end
k = row_of_hour(SS.time);
ok = k >= 1 & k <= n_h;
for name = string(SS.Properties.VariableNames)
    H.(name) = nan(n_h, 1);
    H.(name)(k(ok)) = SS.(name)(ok);
end

LC = read_posted_product(fullfile(products_folder, 'bb2743661p_30_1_lowCloudPeriods.txt'));
H.low_cloud_period = zeros(n_h, 1);
k = row_of_hour(LC.time);
ok = k >= 1 & k <= n_h;
H.low_cloud_period(k(ok)) = LC.LowCloudPeriod(ok);



% -------------------------------------------------------------------
% ------ FM-120 in-cloud droplet number and my s_Nd,FM (hourly) -----
% -------------------------------------------------------------------

in_cloud = fm.lwc_gm3 > opt.FmCloudLwc_gm3 & isfinite(fm.nd_cm3);
k = floor(hours(fm.time - opt.TimeRange(1))) + 1;         % hour [h, h + 1) holds each 5-min start
ok = in_cloud & k >= 1 & k <= n_h;
H.n_fm_cloud = accumarray(k(ok), 1, [n_h 1], @sum, 0);
H.nd_fm_cm3  = accumarray(k(ok), fm.nd_cm3(ok), [n_h 1], @sum, 0) ./ H.n_fm_cloud;   % NaN where 0/0

H.sNdFM_AJB = nan(n_h, 1);
for i = find(isfinite(H.nd_fm_cm3) & all(isfinite(H.ccn_cm3), 2))'
    H.sNdFM_AJB(i) = retrieve_supersaturation_pct(H.nd_fm_cm3(i), H.ccn_cm3(i, :));
end


end
