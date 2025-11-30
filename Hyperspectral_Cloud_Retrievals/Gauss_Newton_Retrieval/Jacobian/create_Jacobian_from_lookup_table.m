%% -- Create Jacobian Matrix from pre computed look up tables --


% By Andrew J. Buggee
%%


function [K,offset] = create_Jacobian_from_lookup_table(inputs,data_inputs,folderName)

% Load the look-up table that is needed for this data set
if strcmp(inputs.data_type,'aviris')==true
    
    error('There are no look-up tables yet!')
    
    
elseif strcmp(inputs.data_type,'modis')==true
    
    model_mean = inputs.model.mean;
    
   
    num_pixels = data_inputs.inputs.pixels.num_2calculate;
    re = data_inputs.inputs.re; % - microns - effective radius - value is in the file name
    tau_c = data_inputs.inputs.tau_c; % - cloud optical depth
    interpGridScaleFactor = data_inputs.inputs.interpGridScaleFactor; % scale factor the will be used to increase the grid size for interpolation.
    
    R = data_inputs.R;
    
    num_bands = size(R,4); % number of bands the look-up table spans
    % lets interpolate the model data to increase our grid
    
    % in our model data, the column space, whixh spans the x direction, varies
    % with tau. The row space, which spans the y direction, vaires with re
    
    [X,Y] = meshgrid(tau_c,re);
    
    % now lets expand the two variables in question to the preferred size
    
    newTau_c = linspace(min(tau_c),max(tau_c),interpGridScaleFactor*length(tau_c));
    new_re = linspace(min(re),max(re),interpGridScaleFactor*length(re));
    
    % set up the new grid to interpolate on
    [Xq,Yq] = meshgrid(newTau_c,new_re);
    
    
    % the jacobian is calculate for each wavelength bin, and for each model
    % parameter
    K = cell(1,num_pixels);
    offset = cell(1,num_pixels);
    for pp = 1:num_pixels
        
        
        % lets shape our data in an easy to use format
        % first extract the data from the bands of interest
        
        % lets compare with non-interpolated data and use the 1km resolution
        % modis reflected data, right now, is in 500 meter resolution. So we
        % have to use the 500 meter pixel indexes
        
        
        R_pix = R(pp,:,:,:);
        R_pix = reshape(R_pix,size(R,2),size(R,3),num_bands);
        
        %offset(:,:,pp)
        % preform 2D interpolation
        newR = zeros(length(new_re),length(newTau_c),num_bands);
%       offset = zeros(1,num_bands);
        for bb = 1:num_bands
            newR(:,:,bb) = interp2(X,Y,R_pix(:,:,bb),Xq,Yq); % new interpolated look up table
%             offset(bb) = interp2(Xq,Yq,newR(:,:,bb),model_mean(2),model_mean(2));
        end
        
        % ----- Find Derivatives with Respect to Each Retrieval Variable -----
        % always taking the first derivatinve here, hence the 1 in diff
        
%        derivative_newR_Y = diff(newR,1,1)./diff(Yq,1,1); % derivative with respect to row space (y -space) which is r
%        derivative_newR_X = diff(newR,1,2)./diff(Xq,1,2); % derivative with respect to column space (x-space) with is tau
        
        % we need these two derivatives to be the same size matrix. So we
        % drop the boundary values for both
        
        % make the derivative matrices the same size

%         derivative_newR_Y(:,end,:) = [];
%         derivative_newR_X(end,:,:) = [];
        
        % average derivatives to fit a plane
        avg_derivative_r = mean(newR(end,:,:) - newR(1,:,:),2)./mean(Yq(end,:) - Yq(1,:)); % derivative in re
        avg_derivative_tau = mean(newR(:,end,:) - newR(:,1,:),1)./mean(Xq(:,end) - Xq(:,1)); % derivative in tau_c
        
        % find the gradient shift graf(F)*x0
        grad_shift  = avg_derivative_r * model_mean(1) + avg_derivative_tau * model_mean(2);
        grad_shift = reshape(grad_shift,num_bands,[]);
        % create (x - x0) 
        newX = repmat(Xq(1:end-1,1:end-1) - model_mean(2),1,1,num_bands); % tau values
        newY = repmat(Yq(1:end-1,1:end-1) - model_mean(1),1,1,num_bands); % re values
        
%         X0 = repmat(Xq(1:end-1,1:end-1),1,1,num_bands);
%         Y0 = repmat(Yq(1:end-1,1:end-1),1,1,num_bands);
        
%         model_by_definition{pp} = reshape(offset,1,1,[]) + derivative_newR_Y.*newY + derivative_newR_X .* newX;
%         model_withOut_shift{pp} = reshape(offset,1,1,[]) + derivative_newR_Y.*Y0 + derivative_newR_X .* X0;
%         model_forced_linear = reshape(offset,1,1,[]) + avg_derivative_r.*newY + avg_derivative_tau .* newX;
         
        offset_try = 0:0.01:1; % values to try for the offset in reflectance to fit the linear plane
        min_val = zeros(1,num_bands);
        index_min = zeros(1,num_bands);
        for bb = 1: num_bands
            sum_squares = zeros(1,length(offset_try));
            for ff = 1:length(offset_try)
            model_forced_linear = offset_try(ff) + avg_derivative_r(bb).*newY(:,:,bb) + avg_derivative_tau(bb) .* newX(:,:,bb);
            sum_squares(ff) = sqrt(sum((newR(1:end-1,1:end-1,bb) - model_forced_linear).^2,'all'));

            end
            
            [min_val(bb),index_min(bb)] = min(sum_squares);
        end
        
        offset{pp} = offset_try(index_min);
        % new forward model = [offset + dF/dr * r + dF/dT * T]
        K{pp} = [(offset{pp}' - grad_shift), reshape(avg_derivative_r,num_bands,1),reshape(avg_derivative_tau,num_bands,1)];

%         for ii = 1:size(derivative_newR_Y,1)
%             for jj = 1:size(derivative_newR_Y,2)
%                 
%                 K{pp}(:,1,ii*jj) = reshape(derivative_newR_Y(ii,jj,:),num_bands,1);
%                 K{pp}(:,2,ii*jj) = reshape(derivative_newR_X(ii,jj,:),num_bands,1);
%                 
%             end
%         end
        
        
        
    end
    
    
end











end