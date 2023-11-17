%% ----- Compare the effectuve radius calculate by modis, and my algorithm -----



% By Andrew J. Buggee
%%

function [] = plot_effRadius_modis_estimates(truth_estimate_table, inputs)


% loop through all rows in inputs.bands2plot

for nn = 1:size(inputs.bands2plot,1)

    if inputs.bands2plot(nn,1)==1 && inputs.bands2plot(nn,2)==7
        % ---------------------------------------------
        % ----------- Plot Bands 1 and 7 --------------
        % ---------------------------------------------

        % extract the modis estimate and my calculation estimates
        modis_R17 = truth_estimate_table.modisR17;

        % The uncertainty is listed as a percent. Lets convert this to a value in
        % microns
        modis_R17_uncert = modis_R17.*(truth_estimate_table.modisR17_uncert./100);          % microns

        % Grab my estimates
        est_R17 = truth_estimate_table.estR17;

        square_diffR17 = truth_estimate_table.squareDiffR17; % the absolute difference between my estimate and the modis estimate
        rms_diff_R17 = sqrt(mean(square_diffR17));

        % find the minimum and maximum values to create a y=x line

        min_est = min(est_R17);
        min_modis = min(modis_R17);

        max_est = max(est_R17);
        max_modis = max(modis_R17);

        min_global = min([min_est,min_modis]);

        max_global = min([max_est,max_modis]);

        x = linspace((0.9 * min_global),(1.1*max_global),150);


        f = figure; plot(x,x,'k-','Linewidth',1)
        hold on; grid on; grid minor
        errorbar(est_R17,modis_R17,modis_R17_uncert,'vertical','m.','MarkerSize',20)
        xlabel('My Estimate: $r_{e}$ $(\mu m)$','Interpreter','latex')
        ylabel('MODIS Estimate: $r_{e}$ $(\mu m)$','Interpreter','latex')
        title(['Bands $1\& 7$ - RMS: ',num2str(rms_diff_R17),' $\mu m$'],'Interpreter','latex')
        set(f, 'Position', [0 0 500 500])
        axis square

        % ---------------------------------------------
        % ----------- Plot Bands 2 and 7 --------------
        % ---------------------------------------------

        % MODIS may have used bands 2 and 7 instead of 1 and 7
        % find the minimum and maximum values to create a y=x line

        % min_est = min(est_R27);
        % min_modis = min(modis_R17);
        %
        % max_est = max(est_R27);
        % max_modis = max(modis_R17);
        %
        % min_global = min([min_est,min_modis]);
        %
        % max_global = min([max_est,max_modis]);
        %
        % x = linspace((0.9 * min_global),(1.1*max_global),150);
        %
        %
        % f = figure; plot(x,x,'k-','Linewidth',1)
        % hold on; grid on; grid minor
        % plot(est_R27,modis_R17,'m.')
        % xlabel('My Estimate: r_{e} (\mum)')
        % ylabel('MODIS Estimate: r_{e} (\mum)')
        % title(['Bands 2&7 - RMS: ',num2str(rms_diff_R27),' \mum'])
        % set(f, 'Position', [0 0 1000 400])



    elseif inputs.bands2plot(nn,1)==1 && inputs.bands2plot(nn,2)==6
        % ---------------------------------------------
        % ----------- Plot Bands 1 and 6 --------------
        % ---------------------------------------------

        modis_R16 = truth_estimate_table.modisR16;
        modis_R16_uncert = modis_R16.*(truth_estimate_table.modisR16_uncert./100);

        est_R16 = truth_estimate_table.estR16;




        square_diffR16 = truth_estimate_table.squareDiffR16; % the absolute difference between my estimate and the modis estimate
        rms_diff_R16 = sqrt(mean(square_diffR16));

        % find the minimum and maximum values to create a y=x line

        min_est = min(est_R16);
        min_modis = min(modis_R16);

        max_est = max(est_R16);
        max_modis = max(modis_R16);

        min_global = min([min_est,min_modis]);

        max_global = min([max_est,max_modis]);

        x = linspace((0.9 * min_global),(1.1*max_global),150);


        f =figure; plot(x,x,'k-','Linewidth',1)
        hold on; grid on; grid minor
        errorbar(est_R16,modis_R16,modis_R16_uncert,'vertical','m.','MarkerSize',10)
        xlabel('My Estimate: r_{e} (\mum)')
        ylabel('MODIS Estimate: r_{e} (\mum)')
        title(['Bands 1&6 - RMS: ',num2str(rms_diff_R16),' \mum'])
        set(f, 'Position', [0 0 1000 400])

    end

end




% ---------------------------------------------
% ----------- Plot Bands 1 and 7 --------------
% ---------------------------------------------

% find the indices of estiamtes that are furthest from their modis
% counterpart

% num2find = 10;
% index_2find = zeros(1,num2find);
%
% for ii = 1:num2find
%
%     [~, index_2find(ii)] = max(square_diffR17);
%
%     square_diffR17(index_2find(ii)) = 0; % set it to a value that will never be chosen!
%
% end




% what do I want to look at with the estimates that deviate the most from
% the modis values?

% first lets plot all pixels again, and then highlight the 10 that deviate
% the most. This allows us to see how they stack up with the full set

% f = figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% plot(est_R17,modis_R17,'m.')
% xlabel('My Estimate: r_{e} (\mum)')
% ylabel('MODIS Estimate: r_{e} (\mum)')
% legend('Perfect Fit','all pixels','Location','best')
% %plot(est_R17(index_2find),modis_R17(index_2find),'c.')
% %legend('Perfect Fit','all pixels',[num2str(num2find),' furthest from line'],'Location','best')
% %title(['Mean Abs Difference: ',num2str(avg_abs_diff),' \mum'])
% title(['Bands 1&7 - RMS: ',num2str(rms_diff_R17),' \mum'])
% set(f, 'Position', [0 0 1000 400])

% find and remove values of tau that modis deems to be greater than 80.
% This is the upper limit I set in my look up tables.
%
% index_80 = truth_estimate_table.modisT17>80;
% est_R17_80 = est_R17(~index_80);
% modis_R17_80 = modis_R17(~index_80);
%
% % redefine the absolute difference since we altered some values up above.
% % Also, ignore any values where the tau is greater than 80
%
% square_diffR17 = truth_estimate_table.squareDiffR17(~index_80); % the absolute difference between my estimate and the modis estimate
% num2find = 10;
% index_2find = zeros(1,num2find);
%
% for ii = 1:num2find
%
%     [~, index_2find(ii)] = max(square_diffR17);
%
%     square_diffR17(index_2find(ii)) = 0; % set it to a value that will never be chosen!
%
% end

% figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% plot(est_R17_80,modis_R17_80,'m.')
% xlabel('My estimate - r_{e} (\mum)')
% ylabel('MODIS estimate - r_{e} (\mum)')
% plot(est_R17_80(index_2find),modis_R17_80(index_2find),'c.')
% legend('Perfect Fit','all pixels',[num2str(num2find),' furthest from line'],'Location','best')
% title(['Mean Abs Difference: ',num2str(mean(abs_diff))])






% -----------------------------------------------------------
% -----------------------------------------------------------
% -----------------------------------------------------------
% PLOT THE SAME GRAPHS BUT FILTER FOR MODIS R_E<24
% THIS IS THE R_E LIMIT IN THE PRE-COMPUTED MIE TABLE
% -----------------------------------------------------------
% -----------------------------------------------------------
% -----------------------------------------------------------
% extract the modis estimate and my calculation estimates

% Find Modis values less than 24 microns, the limit in our pre-computed mie
% table
% index24 = truth_estimate_table.estR17<24;
%
% modis_R17 = truth_estimate_table.modisR17(index24);
% modis_R17_uncert = modis_R17.*(truth_estimate_table.modisR17_uncert(index24)./100);
%
% est_R17 = truth_estimate_table.estR17(index24);
%
% square_diffR17 = truth_estimate_table.squareDiffR17(index24); % the absolute difference between my estimate and the modis estimate
% rms_diff_R17 = sqrt(mean(square_diffR17));
%
% % find the minimum and maximum values to create a y=x line
%
% min_est = min(est_R17);
% min_modis = min(modis_R17);
%
% max_est = max(est_R17);
% max_modis = max(modis_R17);
%
% min_global = min([min_est,min_modis]);
%
% max_global = min([max_est,max_modis]);
%
% x = linspace((0.9 * min_global),(1.1*max_global),150);


% f = figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% errorbar(est_R17,modis_R17,modis_R17_uncert,'vertical','m.','MarkerSize',10)
% xlabel('My Estimate: r_{e} (\mum)')
% ylabel('MODIS Estimate: r_{e} (\mum)')
% title(['Bands 1&7 - RMS: ',num2str(rms_diff_R17),' \mum'])
% set(f, 'Position', [0 0 1000 400])

end