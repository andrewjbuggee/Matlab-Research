function txt = format_value(m, s)
% FORMAT_VALUE  Text for one table cell: "mean ± std", or a single number.
%
%   txt = format_value(m, s)
%
%   Inputs   m   the mean (or a total, or a fraction)
%            s   the standard deviation; NaN when there is none
%   Output   txt  a string, e.g. "0.56 ± 0.27", "323.0 ± 4.8", "473", "0.853",
%                 or "n/a" when m is NaN (no data in that season)
%
%   Rounding rules (the same ones used for every value in the table):
%   * With a standard deviation: the std is given to 2 significant figures,
%     and the mean to the same number of decimal places.
%         0.2714 -> 0.27 (2 decimals)  so 0.5634 -> "0.56 ± 0.27"
%         4.83   -> 4.8  (1 decimal)   so 323.04 -> "323.0 ± 4.8"
%         615.2  -> 615  (0 decimals)  so 902.3  -> "902 ± 615"
%   * Without a standard deviation (hours, fractions): 3 significant figures,
%     trailing zeros removed (472.67 -> "473", 66.17 -> "66.2", 1.000 -> "1").

if ~isfinite(m)
    txt = "n/a";
    return
end

if isfinite(s) && s ~= 0
    % number of decimals that leaves 2 significant figures on the std
    n_dec = max(0, 1 - floor(log10(abs(s))));
    txt = sprintf('%.*f ± %.*f', n_dec, m, n_dec, s);
else
    % 3 significant figures, written without an exponent
    if m == 0
        n_dec = 0;
    else
        n_dec = max(0, 2 - floor(log10(abs(m))));
    end
    txt = sprintf('%.*f', n_dec, m);
    if contains(txt, '.')
        txt = regexprep(txt, '0+$', '');   % drop trailing zeros ...
        txt = regexprep(txt, '\.$', '');   % ... and a dangling decimal point
    end
end
txt = string(txt);

end
