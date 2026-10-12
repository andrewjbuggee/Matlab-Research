function r = make_result(sheet_row, quantity, units, location, product, stats, notes)
% MAKE_RESULT  Package one new line of the table.
%
%   r = make_result(sheet_row, quantity, units, location, product, stats, notes)
%
%   Inputs
%     sheet_row  row number of Kavin's row in the ORIGINAL sheet
%                (Table-EPCAPE_SeasonalAverages.xlsx, sheet Cloud&MetQuantities)
%                that this line checks, e.g. 10 for Liquid Water Path
%     quantity, units, location, product   text for columns A-D of the sheet
%     stats      struct from season_stats.m or season_sums.m (fields mean,
%                std, n; each 1 x 5 = EPCAPE, Spring, Summer, Fall, Winter)
%     notes      text for the Notes column: data, filters and thresholds
%   Output
%     r          struct with all of the above plus r.text, the five cells as
%                they appear in the table (e.g. "0.56 ± 0.27"; see format_value.m)

r.sheet_row = sheet_row;
r.quantity  = string(quantity);
r.units     = string(units);
r.location  = string(location);
r.product   = string(product);
r.mean      = stats.mean;
r.std       = stats.std;
r.n         = stats.n;
r.notes     = string(notes);
r.text      = strings(1, 5);
for k = 1:5
    r.text(k) = format_value(stats.mean(k), stats.std(k));
end

end
