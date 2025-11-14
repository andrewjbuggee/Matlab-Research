function R2 = compute_qqplot_R2(hData, hLine)
% COMPUTE_QQPLOT_R2 Compute R^2 from the handles of a QQ-plot.
%
%   R2 = COMPUTE_QQPLOT_R2(hData) takes hData, the handle to the
%   plotted quantile points (e.g., the first output handle from qqplot),
%   and computes the coefficient of determination between the theoretical
%   quantiles (XData) and the sample quantiles (YData).
%
%   R2 = COMPUTE_QQPLOT_R2(hData, hLine) optionally also takes hLine,
%   the handle of the reference line, though this argument is not used
%   in the computation (included for completeness if you want to inspect
%   the line handle).
%
%   Example usage:
%     h = qqplot(sampleData, 'normal');
%     R2 = compute_qqplot_R2(h(1));
%
%   See also qqplot.

    % Validate input
    if nargin < 1
        error('You must supply at least one handle (the data points) as input.');
    end
    
    % Extract the quantile points handle
    % We assume hData is the handle whose XData are the theoretical quantiles
    % and YData are the sample quantiles.
    % If hData is an array of handles, we pick the first one.
    if numel(hData) > 1
        hPts = hData(1);
    else
        hPts = hData;
    end
    
    % Get x and y values
    x = get(hPts, 'XData');
    y = get(hPts, 'YData');
    
    % Remove any NaN or Inf if present
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    
    % If too few points, warn
    if numel(x) < 2
        warning('Less than 2 valid quantile‐points found; R^2 is set to NaN.');
        R2 = NaN;
        return;
    end
    
    % Fit a linear regression y = a + b*x
    coeffs = polyfit(x, y, 1);
    yFit = polyval(coeffs, x);
    
    % According to MathWorks docs: 
    %   SST = total sum of squares = sum((y - mean(y))^2)
    %   SSE = sum of squared errors (residuals) = sum((y - yFit)^2)
    %   R^2 = 1 – SSE/SST.  :contentReference[oaicite:3]{index=3}
    SSE = sum((y - yFit).^2);
    SST = sum((y - mean(y)).^2);
    
    % Protect against zero‐variance y
    if SST == 0
        warning('Zero variance in Y-data: SST=0; R^2 is set to NaN.');
        R2 = NaN;
    else
        R2 = 1 - SSE/SST;
    end
    
    % % (Optional) annotate the plot
    % try
    %     ax = ancestor(hPts, 'axes');
    %     xlimVals = xlim(ax);
    %     ylimVals = ylim(ax);
    %     % Choose a location for text: e.g., 5% in from left and 10% down from top
    %     xLoc = xlimVals(1) + 0.05 * diff(xlimVals);
    %     yLoc = ylimVals(2) - 0.10 * diff(ylimVals);
    %     text(ax, xLoc, yLoc, sprintf('R^2 = %.4f', R2), ...
    %          'FontWeight','bold','Color','b', 'VerticalAlignment','top');
    % catch
    %     % If annotation fails, ignore
    % end
    % 
    % % Display result
    % fprintf('Computed R^2 (coefficient of determination) = %.4f\n', R2);
end


