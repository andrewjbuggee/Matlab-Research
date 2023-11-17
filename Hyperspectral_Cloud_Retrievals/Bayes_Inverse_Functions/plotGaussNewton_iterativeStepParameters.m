%% Plot Gauss-Newton Iterative Solution Parameters to Learn how the solution changes with each iteration


% By Andrew J. Buggee
%%

function f = plotGaussNewton_iterativeStepParameters(GN_inputs, GN_outputs)

num_iterations = GN_inputs.GN_iterations;
num_pixels = GN_inputs.numPixels2Calculate;
num_parameters = GN_inputs.num_model_parameters;


% Define left and right color axes
leftAxColor = [0, 0, 0];    % black
rightAxColor = [0, 0, 1];   % blue

if num_pixels <=3
    
    for pp = 1:num_pixels
        
        f(pp) = figure;
        
        subplot(1,2,1)
        set(f(pp), 'defaultAxesColorOrder', [leftAxColor; rightAxColor]);
        for xx = 1:num_parameters
            if xx==1
                yyaxis left
                plot((1:num_iterations+1), GN_outputs.retrieval(1,:,pp),'k.')
                hold on
                % Plot the a priori value for r_top, which is also that of
                % r_bottom
                yline(GN_inputs.model.mean(1),'k',{'r_{e} a priori'})
                
            elseif xx==2
                plot((1:num_iterations+1), GN_outputs.retrieval(2,:,pp),'k*', 'MarkerSize',3)
            else
                yyaxis right
                plot((1:num_iterations+1), GN_outputs.retrieval(3,:,pp),'b.')
                ylabel('Optical Depth')
                hold on
                % Plot the a priori value as a constant line
                yline(GN_inputs.model.mean(3), 'b', {'\tau_c a priori'});
            end
        end
        
        yyaxis left
        
        
        ylim([min(GN_outputs.retrieval(:,2))*0.9, max(GN_outputs.retrieval(:,1))*1.1])
        grid on; grid minor
        xlabel('Iteration')
        ylabel('Retrieval Value')
        set(f(pp), 'Position', [0 0 1000 400])
        
        
    end
    
    
    
elseif num_pixels > 3
    
    % if there are a bunch of pixels, we will just grab a random subset of
    % 3 to plot
    rand_index = randsample(num_pixels,3); % random sampling without replacement
    
    for pp = 1:length(rand_index)
        
        
    end
    
    
    
end





end





