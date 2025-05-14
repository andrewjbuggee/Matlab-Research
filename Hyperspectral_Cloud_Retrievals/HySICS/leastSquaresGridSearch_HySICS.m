%% ----- Reflectance Function - Coarse Grid Search -----


% right now, this can only use two wavelengths to perform a minimum least
% squares grid search

% By Andrew J. Buggee

%%

function [tblut_retrieval] = leastSquaresGridSearch_HySICS(simulated_reflectance, modelRefl, inputs)



% extract inputs
re = inputs.RT.re;
tau_c = inputs.RT.tau_c;
interpGridScaleFactor = inputs.interpGridScaleFactor;
bands2run = inputs.bands2run_from_set_of_measurements; % bands to run through uvspec

% extract size of forward modeled calcualtions
num_rEff = length(re);
num_tauC = length(tau_c);
num_wl = 2;



% Define the midpoint of each spectral channel
band_midPoint = mean(inputs.RT.wavelengths2run, 2);      % nm


% lets interpolate the model data to increase our grid

% in our model data, the column space, which spans the x direction, varies
% with tau. The row space, which spans the y direction, vaires with re

[Re0, T0] = meshgrid(re,tau_c);

% now lets expand the two variables in question to the preferred size

newTau_c = linspace(min(tau_c),max(tau_c),interpGridScaleFactor*length(tau_c));
new_re = linspace(min(re),max(re),interpGridScaleFactor*length(re));




% set up the new grid to interpolate on
[Re, T] = meshgrid(new_re, newTau_c);



% reshape the libRadTran forward modeled reflectance to be an array where one
% dimension varies with optical depth and another varies with
% effective radius
modelRefl_band1_array = reshape(modelRefl(1:2:length(modelRefl)), num_tauC, num_rEff);
modelRefl_band2_array = reshape(modelRefl(2:2:length(modelRefl)), num_tauC, num_rEff);



% preform 2D interpolation
interp_modelRefl_band1 = interp2(Re0, T0, modelRefl_band1_array, Re, T);
interp_modelRefl_band2 = interp2(Re0, T0, modelRefl_band2_array, Re, T);




% finally, lets create an array with the same shape as the new
% interpolated array, but where every value is the observed
% reflectance
% for band 1...
observations_newGrid_band1 = repmat(simulated_reflectance(bands2run(1)),...
                  length(newTau_c), length(new_re));
% for band 2...
observations_newGrid_band2 = repmat(simulated_reflectance(bands2run(2)),...
                  length(newTau_c), length(new_re));

%% ---- lets view the surfaces of the model -----

band2Plot = 2;

if inputs.flags.plotMLS_figures == true
    surfPlots4modisModel_andObs(T,Re,interp_modelRefl(:,:,band2Plot),observations_newGrid(:,:,band2Plot),band_midPoint(1,band2Plot))
end
%% ----- Least Squares Difference ------

% take the difference between the model and the observaions at each
% wavelength. The we square the difference and sum along wavelenghts. Then
% we take the sqrt

leastSquaresGrid = sqrt(mean((cat(3, interp_modelRefl_band1, interp_modelRefl_band2) - ...
    cat(3, observations_newGrid_band1, observations_newGrid_band2)).^2,3));

% If we assume that there is a unique global minimum, we can search for it
% using the minimum function. If we have more complicated scenario we need
% to use the optimization tools to search for global and local minima

[tblut_retrieval.minLSD,index] = min(leastSquaresGrid,[],'all','linear');
[row,col] = ind2sub(size(leastSquaresGrid),index);

% Save the effective radius and optical depth associated with the
% minimum RMS difference
tblut_retrieval.minRe = Re(row,col);
tblut_retrieval.minTau = T(row,col);

% Save the reflectance associated with the minimum RMS difference
tblut_retrieval.reflectance = [interp_modelRefl_band1(row, col); interp_modelRefl_band2(row, col)];


% lets look at the least squares grid
if inputs.flags.plotMLS_figures == true

    figure; s = surf(T,Re,leastSquaresGrid);
    xlabel('Optical Depth')
    ylabel('Effective Radius (\mum)')
    zlabel(['Least Squares Difference'])
    s.EdgeColor = 'k';
    s.EdgeAlpha = 0.1;
    colorbar
    set(gcf,"Position", [0 0 1000 700])
    title('Error Function')


    % -- Create contour plot of RMS residual ---

    % Create figure
    figure;
    colormap(hot);

    % Create axes
    axes1 = axes;
    hold(axes1,'on');

    % plot the whole space?
    T_idx = 250;
    % number of levels to plot
    n = 15;

    % Create contour
    [c1,h1] = contour(T(:, 1:T_idx), Re(:, 1:T_idx), leastSquaresGrid(:, 1:T_idx), n,'LineWidth',3);
    clabel(c1,h1,'FontSize',20,'FontWeight','bold');

    % round off the level list
    h1.LevelList = round(h1.LevelList,4);  %rounds levels to 3rd decimal place
    clabel(c1,h1)

    % plot the global minimum
    hold on; plot(tblut_retrieval.minTau(pp), tblut_retrieval.minRe(pp), 'k.', 'MarkerSize', 25)

    % make a legend
    legend('RMS differences', 'Global Minimum', 'Location', 'best', 'Interpreter', 'latex');

    % Create ylabel
    xlabel('$\tau_{c}$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    % Create xlabel
    ylabel('$r_{e}$ $(\mu m)$','FontWeight','bold','Interpreter','latex', 'Fontsize', 35);

    % Create title
    title('RMS difference between EMIT and libRadTran calculations','Interpreter','latex');

    box(axes1,'on');
    grid(axes1,'on');
    axis(axes1,'tight');
    hold(axes1,'off');
    % Set the remaining axes properties
    set(axes1,'BoxStyle','full','Layer','top','XMinorGrid','on','YMinorGrid','on','ZMinorGrid',...
        'on');
    % Create colorbar
    colorbar(axes1);
    % set the figure size to be proportional to the length of the r_top and
    % r_bot vectors
    %set(gcf, 'Position', [0 0 1200, 1200*(length(r_bot)/length(r_top))])
    set(gcf, 'Position', [0 0 900 900])



end




save(inputs.save_mat_filename,"tblut_retrieval",'-append'); % save inputSettings to the same folder as the input and output file



end









