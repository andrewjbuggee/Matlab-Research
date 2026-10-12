function results = calc_rows27_30_activated_fraction(kavin_folder)
% CALC_ROWS27_30_ACTIVATED_FRACTION  Rows 27 and 30: cloud and haze activated fraction.
%
%   results = calc_rows27_30_activated_fraction(kavin_folder)
%
%   DATA: ET_Dcrit_nativetime.csv (from Abbey Williams' code, via Kavin), one
%   row per sample, ~1-4 s apart:
%     time             e.g. "10-Mar-2023 22:38:16" (UTC)
%     Dta              a droplet-residual diameter
%     Dta_critical     the critical diameter for activation
%     ClrCldHaz_class  1 = cloud, 2 = haze, empty = neither
%   (Neither diameter is defined in the file; they come from Abbey's code.)
%
%   WHAT IS COMPUTED (Kavin's count, applied in each season):
%     take the rows of one class that have BOTH diameters; the activated
%     fraction = share of those rows with Dta > Dta_critical. Every row counts
%     once (samples, not time). No ± is given for a fraction.
%
%   LIMITED DATA: Dta_critical exists only in Jun-Sep and Dec 2023, so
%   "Spring" is June alone, "Fall" is December alone, and there is no
%   winter value (and no fall or winter haze value).

path = fullfile(kavin_folder, 'ET_Dcrit_nativetime.csv');
opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'time', 'datetime');
opts = setvaropts(opts, 'time', 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'DatetimeLocale', 'en_US');
opts = setvartype(opts, {'Dta', 'Dta_critical', 'ClrCldHaz_class'}, 'double');
T = readtable(path, opts);

masks = sheet_season_masks(T.time);
names = ["Cloud Droplet Activated Fraction", "Haze Droplet Activated Fraction"];
rows  = [27, 30];
results = [];
for code = 1:2
    use = (T.ClrCldHaz_class == code) & isfinite(T.Dta) & isfinite(T.Dta_critical);
    activated = double(T.Dta > T.Dta_critical);   % 1 = activated, 0 = not
    activated(~use) = NaN;                        % rows not used
    stats = season_stats(activated, masks);       % the mean of 0/1 values = the fraction
    stats.std(:) = NaN;                           % no ± for a fraction
    results = [results, make_result(rows(code), names(code), "fraction", "Soledad", ...
        "ET_Dcrit_nativetime.csv (Abbey Williams)", stats, ...
        sprintf("Kavin's count re-run on his file, per season. filters: ClrCldHaz_class = %d; rows with both ", code) + ...
        "Dta and Dta_critical; fraction with Dta > Dta_critical (per sample, not time-weighted); sheet seasons " + ...
        "(UTC). Not independent: Dta and Dta_critical come from Abbey's code. Dta_critical only Jun-Sep and Dec 2023.")]; %#ok<AGROW>
end

end
