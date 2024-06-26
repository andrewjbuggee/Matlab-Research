%% ----- Reflectance Function - Coarse Grid Search -----


% right now, this can only use two wavelengths to perform a minimum least
% squares grid search

% By Andrew J. Buggee

%%

function [minVals] = leastSquaresGridSearch_EMIT(emitRefl, modelRefl,inputs)



% extract inputs
re = inputs.RT.re;
tau_c = inputs.RT.tau_c;
interpGridScaleFactor = inputs.interpGridScaleFactor;
bands2run = inputs.bands2run; % bands to run through uvspec

% folder where the .mat file of reflectances is stored
saveCalcs_folder = inputs.folder2save.reflectance_calcs;
% save calculations
saveCalcs_filename = inputs.reflectance_calculations_fileName;


numPixels = size(modelRefl,1); % number of pixels to preform grid search on



% Define the midpoint of each spectral channel
band_midPoint = inputs.source.wavelength(inputs.bands2run);      % nm


% lets interpolate the model data to increase our grid

% in our model data, the column space, which spans the x direction, varies
% with tau. The row space, which spans the y direction, vaires with re

[T0,Re0] = meshgrid(tau_c,re);

% now lets expand the two variables in question to the preferred size

newTau_c = linspace(min(tau_c),max(tau_c),interpGridScaleFactor*length(tau_c));
new_re = linspace(min(re),max(re),interpGridScaleFactor*length(re));

% set up the new grid to interpolate on
[T,Re] = meshgrid(newTau_c,new_re);



for pp = 1:numPixels


    % lets shape our data in an easy to use format

    % grab the observations for the pair of bands desired for the
    % retrieval
    observations = emitRefl.value(bands2run', :);     % 1/sr - reflectance

    if iscell(modelRefl) == true

        modelRefl_array = [modelRefl{pp,:,:, :}];
        modelRefl_array = reshape(modelRefl_array,size(modelRefl,2),size(modelRefl,3));



    elseif isnumeric(modelRefl) == true

        modelRefl_array = modelRefl(pp,:,:,:);

        % reshape the libRadTran reflectance to be an array where one
        % dimension varies with optical depth and another varies with
        % effective radius
        modelRefl_array = reshape(modelRefl_array, size(modelRefl,2), size(modelRefl,3), length(bands2run));


    end

    % preform 2D interpolation
    interp_modelRefl = zeros(length(new_re),length(newTau_c),size(modelRefl_array,3));

    % the third dimension is the spectral dimension. Step through each
    % wavelength and perform the 2D interpolation
    for bb = 1:length(bands2run)
        interp_modelRefl(:,:,bb) = interp2(T0,Re0,modelRefl_array(:,:,bb),T,Re);
    end

    % finally, lets create an array with the same shape as the new
    % interpolated array, but where every value is the observed
    % reflectance
    % for band 1...
    observations_newGrid(:,:,1) = repmat(observations(1,pp),length(new_re),length(newTau_c));
    % for band 2...
    observations_newGrid(:,:,2) = repmat(observations(2, pp),length(new_re),length(newTau_c));

    %% ---- lets view the surfaces of the model -----

    band2Plot = 2;

    if inputs.flags.plotMLS_figures == true
        surfPlots4modisModel_andObs(T,Re,interp_modelRefl(:,:,band2Plot),observations_newGrid(:,:,band2Plot),band_midPoint(1,band2Plot))
    end
    %% ----- Least Squares Difference ------

    % take the difference between the model and the observaions at each
    % wavelength. The we square the difference and sum along wavelenghts. Then
    % we take the sqrt

    leastSquaresGrid = sqrt(mean((interp_modelRefl - observations_newGrid).^2,3));

    % If we assume that there is a unique global minimum, we can search for it
    % using the minimum function. If we have more complicated scenario we need
    % to use the optimization tools to search for global and local minima

    [minVals.minLSD(pp),index] = min(leastSquaresGrid,[],'all','linear');
    [row,col] = ind2sub(size(leastSquaresGrid),index);
    
    % Save the effective radius and optical depth associated with the
    % minimum RMS difference
    minVals.minRe(pp) = Re(row,col);
    minVals.minTau(pp) = T(row,col);

    % Save the reflectance associated with the minimum RMS difference
    minVals.reflectance(:, pp) = [interp_modelRefl(row, col, 1); interp_modelRefl(row, col, 2)];


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
        hold on; plot(minVals.minTau(pp), minVals.minRe(pp), 'k.', 'MarkerSize', 25)

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




    save([saveCalcs_folder, saveCalcs_filename],"minVals",'-append'); % save inputSettings to the same folder as the input and output file



end


end






