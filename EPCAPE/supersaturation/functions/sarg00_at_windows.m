function s_pct = sarg00_at_windows(modes30, w_ms, T_K, p_Pa, kappa_sub, kappa_coarse, formula)
% SARG00_AT_WINDOWS  Parcel supersaturation s_ARG00 (%) for a list of 30-min windows.
%
%   s_pct = sarg00_at_windows(modes30, w_ms, T_K, p_Pa, kappa_sub, kappa_coarse, formula)
%
%   CREDIT: method of Jeramy L. Dedrick (Dedrick et al. 2026, J. Atmos. Sci., Sec. 2d,
%   p. 1066, supplement Text S10 and Table S7). His driver code for this step is not
%   public; this is a reconstruction from the paper.
%
%   For each window, builds the three Pier lognormal modes (coarse, accumulation,
%   Aitken; the row order of smax_arg2000_pct and of his ParcelActivation_ARG2000.m)
%   and evaluates the Abdul-Razzak & Ghan (2000) maximum supersaturation with the
%   window's free convective velocity w*, cloud-base T and p. Per supplement
%   Text S10 (p. 8), the Aitken and accumulation modes share the submicron Pier
%   kappa and the coarse mode takes kappa = 0.87.
%
%   Inputs (one row per window)
%     modes30       table with the Pier mode fits averaged over the window:
%                   Pier_{ait,acc,coarse}_N (cm-3), Pier_{ait,acc,coarse}_diam (um, median
%                   diameter), Pier_{ait,acc,coarse}_gsd (geometric standard deviation)
%     w_ms          updraft (w*), m s-1
%     T_K, p_Pa     cloud-base temperature (K) and pressure (Pa)
%     kappa_sub     kappa of the Aitken and accumulation modes (scalar or one per window)
%     kappa_coarse  kappa of the coarse mode (scalar)
%     formula       "textbook" (smax_arg2000_pct, ARG2000 eq. 6 as published) or
%                   "posted" (smax_posted_arg2000_pct, a line-by-line port of
%                   Jeramy's posted ParcelActivation_ARG2000.m, issues included)
%   Output
%     s_pct         S_max, % (NaN where an input is missing or w <= 0)

n = numel(w_ms);
if isscalar(kappa_sub), kappa_sub = repmat(kappa_sub, n, 1); end
switch formula
    case "textbook", smax = @smax_arg2000_pct;
    case "posted",   smax = @smax_posted_arg2000_pct;
    otherwise, error('formula must be "textbook" or "posted", not "%s"', formula);
end

s_pct = NaN(n, 1);
for i = 1:n
    % rows: coarse, accumulation, Aitken; columns: N (m-3), D_g (m), sigma_g, kappa
    modes = [modes30.Pier_coarse_N(i) * 1e6, modes30.Pier_coarse_diam(i) * 1e-6, modes30.Pier_coarse_gsd(i), kappa_coarse;
             modes30.Pier_acc_N(i) * 1e6,    modes30.Pier_acc_diam(i) * 1e-6,    modes30.Pier_acc_gsd(i),    kappa_sub(i);
             modes30.Pier_ait_N(i) * 1e6,    modes30.Pier_ait_diam(i) * 1e-6,    modes30.Pier_ait_gsd(i),    kappa_sub(i)];
    if any(~isfinite(modes(:))) || ~isfinite(w_ms(i)) || w_ms(i) <= 0 || ~isfinite(T_K(i)) || ~isfinite(p_Pa(i))
        continue
    end
    s_pct(i) = smax(modes, w_ms(i), T_K(i), p_Pa(i));
end

end
