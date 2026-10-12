function results = calc_rows_15min_file(F)
% CALC_ROWS_15MIN_FILE  Rows computed with Kavin's method on Abbey Williams' 15-min file.
%
%   results = calc_rows_15min_file(F)        F = read_15min_file(...)
%
%   Kavin's AMSforAndrew.mlx works on Abbey Williams' 15-min file,
%   EPCAPE_15mindata.nc. Its ClrCldHaz_class is 1 (cloud) or 2 (haze), set
%   only while the CVI inlet was on. His method, re-run here unrounded:
%     * a quantity's value = mean ± std over the 15-min intervals of the class;
%     * hours = number of 15-min intervals of the class / 4;
%     * seasons = calendar quarters of the UTC time (winter = Jan-Mar);
%       EPCAPE = every interval in the file.
%   Sheet rows:
%     25 / 28  cloud / haze droplet number concentration  FM120_Nd (cm-3)
%     26 / 29  cloud / haze droplet effective diameter    sum(N D^3) / sum(N D^2)
%              of FM120_DSD, with each bin's UPPER edge as its diameter D
%              (that is what reproduces Kavin's values)
%     33 / 34  cloud / haze GCVI-AMS sampling hours (two lines: Kavin's count,
%              and only the intervals in which the AMS reported, sheet seasons)
%     35-42    cloud (35-38) and haze (39-42) residual sulfate, nitrate,
%              ammonium, organics: AMS_mass (ug m-3; already divided by the
%              GCVI enhancement factor). These are concentrations, not mass
%              fractions as the sheet's labels say.
%     54 / 55  cloud / haze liquid water content          FM120_LWC (g m-3)

q_masks = calendar_quarter_masks(F.time);      % Kavin's seasons (UTC months)
is_cloud = (F.class == 1);
is_haze  = (F.class == 2);
src = "EPCAPE_15mindata.nc (Abbey's 15-min file)";
kavin_note = "Kavin's method re-run on his file (unrounded).";

% Effective diameter of each 15-min size distribution, bins at their UPPER edges:
%   D_eff = sum_i N_i D_i^3 / sum_i N_i D_i^2      (um)
D = F.diam_um(:).';                            % 1 x 30 bin diameters
d_eff_um = sum(F.dsd .* D.^3, 2) ./ sum(F.dsd .* D.^2, 2);
d_eff_um(~isfinite(d_eff_um)) = NaN;           % no droplets -> no diameter

% Helper: Kavin's mean ± std of `values` over the intervals selected by `sel`
kavin_stats = @(values, sel) season_stats(values .* nan_unless(sel), q_masks);

results = [];
% --- Rows 25, 26, 28, 29: droplet number and effective diameter --------------------------
results = [results, make_result(25, "Cloud Droplet Number Concentration", "cm-3", "Soledad", src, ...
    kavin_stats(F.nd_cm3, is_cloud), kavin_note)];
results = [results, make_result(26, "Cloud Droplet Effective Diameter", "um", "Soledad", src, ...
    kavin_stats(d_eff_um, is_cloud), kavin_note)];
results = [results, make_result(28, "Haze Droplet Number Concentration", "cm-3", "Soledad", src, ...
    kavin_stats(F.nd_cm3, is_haze), kavin_note)];
results = [results, make_result(29, "Haze Droplet Effective Diameter", "um", "Soledad", src, ...
    kavin_stats(d_eff_um, is_haze), kavin_note)];

% --- Rows 33, 34: sampling hours ----------------------------------------------------------
% the AMS reported in an interval if any of its 5 species has a value
has_ams = any(isfinite(F.ams_ugm3), 2);
s_masks = sheet_season_masks(F.time);
names = ["On-Site Cloud Sampling for AMS Residuals", "On-Site Haze Sampling for AMS Residuals"];
for k = 1:2
    sel = (F.class == k);
    results = [results, make_result(32 + k, names(k), "hr", "Soledad", src, ...
        season_sums(0.25 * sel, q_masks), "Kavin's code re-run on his file, not rounded.")]; %#ok<AGROW>
    results = [results, make_result(32 + k, names(k), "hr", "Soledad", "EPCAPE_15mindata.nc classes + AMS", ...
        season_sums(0.25 * (sel & has_ams), s_masks), ...
        sprintf("filters: Abbey's 15-min intervals with ClrCldHaz_class = %d in which the AMS reported ", k) + ...
        "(AMS_mass present), x 15 min; sheet seasons (UTC). Kavin's count also includes the class intervals " + ...
        "without AMS data.")]; %#ok<AGROW>
end

% --- Rows 35-42: residual concentrations ------------------------------------------------------
sheet_species = ["sulfate", "nitrate", "ammonium", "organics"];      % order of rows 35-38 (and 39-42)
label_species = ["Sulfate", "Nitrate", "Ammonium", "Organic"];
for k = 1:2                                                            % 1 = cloud, 2 = haze
    for s = 1:4
        col = find(F.species == sheet_species(s));                     % column of AMS_mass
        row = 34 + 4 * (k - 1) + s;                                    % 35..38, then 39..42
        if k == 1
            quantity = "Cloud Droplet Residual NR " + label_species(s) + " Mass Fraction";
        else
            quantity = "Haze Droplet Residual NR " + label_species(s) + " Mass Fraction";
        end
        results = [results, make_result(row, quantity, "ug m-3", "Soledad", src, ...
            kavin_stats(F.ams_ugm3(:, col), F.class == k), "Kavin's code re-run on his file.")]; %#ok<AGROW>
    end
end

% --- Rows 54, 55: liquid water content ------------------------------------------------------
results = [results, make_result(54, "Cloud Liquid Water Content", "g m-3", "Soledad", src, ...
    kavin_stats(F.lwc_gm3, is_cloud), kavin_note)];
results = [results, make_result(55, "Haze Liquid Water Content", "g m-3", "Soledad", src, ...
    kavin_stats(F.lwc_gm3, is_haze), kavin_note)];

end
