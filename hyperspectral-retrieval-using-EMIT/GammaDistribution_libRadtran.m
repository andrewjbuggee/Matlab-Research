classdef GammaDistribution_libRadtran < prob.ToolboxFittableParametricDistribution
% GammaDistribution_libRadtran - Custom gamma distribution following libRadTran definition
%
%    An object of the GammaDistribution_libRadtran class represents a gamma
%    probability distribution with effective radius (r_eff) and shape 
%    parameter (alpha) as defined in libRadTran manual equation 6.10.
%
%    GammaDistribution_libRadtran methods:
%       cdf                   - Cumulative distribution function
%       fit                   - Fit distribution to data
%       icdf                  - Inverse cumulative distribution function
%       mean                  - Mean
%       median                - Median
%       pdf                   - Probability density function
%       random                - Random number generation
%       std                   - Standard deviation
%       var                   - Variance
%
%    GammaDistribution_libRadtran properties:    
%       r_eff                 - Effective radius parameter
%       alpha                 - Shape parameter

    properties(Constant)
        DistributionName = 'GammaLibRadtran';
    end

    properties(Dependent=true)
        r_eff  % Effective radius
        alpha  % Shape parameter
    end
    
    properties(Constant)
        NumParameters = 2;
        ParameterNames = {'r_eff' 'alpha'};
        ParameterDescription = {'effective radius' 'shape parameter'};
    end

    properties(GetAccess='public',SetAccess='protected')
        ParameterValues
    end

    methods
        function pd = GammaDistribution_libRadtran(r_eff, alpha)
            if nargin==0
                r_eff = 10;  % Default effective radius
                alpha = 7;   % Default shape parameter
            end
            checkargs(r_eff, alpha);
            
            pd.ParameterValues = [r_eff alpha];
            pd.ParameterIsFixed = [true true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
        
        function m = mean(this)
            % Mean of the gamma distribution
            r_eff = this.r_eff;
            alpha = this.alpha;
            b = gamma(4 + alpha) / (r_eff * gamma(3 + alpha));
            m = (alpha + 1) / b;
        end
        
        function v = var(this)
            % Variance of the gamma distribution
            r_eff = this.r_eff;
            alpha = this.alpha;
            b = gamma(4 + alpha) / (r_eff * gamma(3 + alpha));
            v = (alpha + 1) * (alpha + 2) / (b^2) - this.mean()^2;
        end
        
        function s = std(this)
            s = sqrt(this.var());
        end
    end
    
    methods
        function this = set.r_eff(this, r_eff)
            checkargs(r_eff, this.alpha);
            this.ParameterValues(1) = r_eff;
            this = invalidateFit(this);
        end
        
        function this = set.alpha(this, alpha)
            checkargs(this.r_eff, alpha);
            this.ParameterValues(2) = alpha;
            this = invalidateFit(this);
        end
        
        function r_eff = get.r_eff(this)
            r_eff = this.ParameterValues(1);
        end
        
        function alpha = get.alpha(this)
            alpha = this.ParameterValues(2);
        end
    end
    
    methods(Static)
        function pd = fit(x, varargin)
            [x, cens, freq] = prob.ToolboxFittableParametricDistribution.processFitArgs(x, varargin{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(x, cens, freq, 'GammaLibRadtran');
            freq = ones(size(x));

            % Estimate parameters using method of moments
            % Mean = (alpha + 1) / b, where b = gamma(4+alpha)/(r_eff*gamma(3+alpha))
            % This requires iterative solution; use simple initial estimates
            m = mean(x);
            v = var(x);
            
            % Initial guess: use standard gamma distribution relationships
            alpha_init = max(0.1, (m^2 / v) - 1);
            r_eff_init = m;
            
            % Refine using maximum likelihood
            options = statset('MaxIter', 1000, 'TolX', 1e-6);
            params0 = [r_eff_init, alpha_init];
            
            objfun = @(params) prob.GammaDistribution_libRadtran.likefunc(params, x);
            params = fminsearch(objfun, params0, options);
            
            r_eff = max(0.001, params(1));
            alpha = max(0.001, params(2));

            pd = prob.GammaDistribution_libRadtran(r_eff, alpha);
            pd.ParameterIsFixed = [false false];
            
            [nll, acov] = prob.GammaDistribution_libRadtran.likefunc([r_eff alpha], x);
            pd.ParameterCovariance = acov;
            pd.NegativeLogLikelihood = nll;
            pd.InputData = struct('data', x, 'cens', [], 'freq', freq);
        end

        function [nll, acov] = likefunc(params, x)
            r_eff = params(1);
            alpha = params(2);
            
            if r_eff <= 0 || alpha < 0
                nll = Inf;
                acov = zeros(2);
                return;
            end
            
            pdf_vals = prob.GammaDistribution_libRadtran.pdffunc(x, r_eff, alpha);
            pdf_vals(pdf_vals <= 0) = realmin;
            
            nll = -sum(log(pdf_vals));
            
            % Approximate covariance using finite differences
            n = length(x);
            acov = eye(2) * (1/n);
        end
        
        function y = pdffunc(x, r_eff, alpha)
            % PDF following libRadTran definition
            b = gamma(4 + alpha) / (r_eff * gamma(3 + alpha));
            N = 1 / (b^(-1 - alpha) * gamma(1 + alpha));
            
            y = N * x.^alpha .* exp(-b .* x);
            y(x < 0) = 0;
            y(~isfinite(y)) = 0;
        end
        
        function y = cdffunc(x, r_eff, alpha)
            % CDF using numerical integration
            b = gamma(4 + alpha) / (r_eff * gamma(3 + alpha));
            
            y = zeros(size(x));
            for i = 1:numel(x)
                if x(i) <= 0
                    y(i) = 0;
                else
                    % Use incomplete gamma function
                    y(i) = gammainc(b * x(i), alpha + 1);
                end
            end
            y(isnan(x)) = NaN;
        end
        
        function y = invfunc(p, r_eff, alpha)
            % Inverse CDF using numerical solver
            if nargin < 2, r_eff = 10; end
            if nargin < 3, alpha = 7; end
            
            y = zeros(size(p));
            for i = 1:numel(p)
                if p(i) < 0 || p(i) > 1
                    y(i) = NaN;
                elseif p(i) == 0
                    y(i) = 0;
                elseif p(i) == 1
                    y(i) = Inf;
                else
                    % Use fzero to find inverse
                    fun = @(x) prob.GammaDistribution_libRadtran.cdffunc(x, r_eff, alpha) - p(i);
                    y(i) = fzero(fun, [0, 10*r_eff]);
                end
            end
        end
        
        function y = randfunc(r_eff, alpha, varargin)
            % Generate random numbers
            u = rand(varargin{:});
            y = prob.GammaDistribution_libRadtran.invfunc(u, r_eff, alpha);
        end
    end
    
    methods(Static, Hidden)
        function info = getInfo
            info = getInfo@prob.ToolboxFittableParametricDistribution('prob.GammaDistribution_libRadtran');
            
            info.name = prob.GammaDistribution_libRadtran.DistributionName;
            info.code = 'gammalibradtran';
            info.support = [0, Inf];
            info.closedbound = [true false];
            info.iscontinuous = true;
            info.islocscale = false;
            info.uselogpp = false;
            info.optimopts = true;
            info.logci = [true true];
        end
    end
end

function checkargs(r_eff, alpha)
if ~(isscalar(r_eff) && isnumeric(r_eff) && isreal(r_eff) && r_eff > 0 && isfinite(r_eff))
    error('r_eff must be a positive finite numeric scalar.')
end
if ~(isscalar(alpha) && isnumeric(alpha) && isreal(alpha) && alpha >= 0 && isfinite(alpha))
    error('alpha must be a non-negative finite numeric scalar.')
end
end