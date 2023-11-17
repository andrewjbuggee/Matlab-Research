%% ----- Compare the cloud optical depth calculate by modis, and my algorithm -----


% By Andrew J. Buggee
%%

function [] = plot_tau_modis_estimates(truth_estimate_table, inputs)

% loop through all rows in inputs.bands2plot

for nn = 1:size(inputs.bands2plot,1)

    if inputs.bands2plot(nn,1)==1 && inputs.bands2plot(nn,2)==7




        % ---------------------------------------------
        % ----------- Plot Bands 1 and 7 --------------
        % ---------------------------------------------

        % extract the modis estimate and my calculation estimates
        modis_T17 = truth_estimate_table.modisT17;

        % The uncertainty is listed as a percentage. Lets convert it to a
        % reflectance
        modis_T17_uncert = modis_T17.*(truth_estimate_table.modisT17_uncert./100);
        est_T17 = truth_estimate_table.estT17;
        square_diffT17 = truth_estimate_table.squareDiffT17; % the absolute difference between my estimate and the modis estimate

        rms_diff_T17 = sqrt(mean(square_diffT17));


        % find the minimum and maximum values to create a y=x line

        min_est = min(est_T17);
        min_modis = min(modis_T17);

        max_est = max(est_T17);
        max_modis = max(modis_T17);

        min_global = min([min_est,min_modis]);

        max_global = min([max_est,max_modis]);

        x = linspace((0.9 * min_global),(1.1*max_global),150);


        f = figure; plot(x,x,'k-','Linewidth',1)
        hold on; grid on; grid minor
        errorbar(est_T17,modis_T17, modis_T17_uncert,'vertical','m.','MarkerSize',20)
        xlabel('My Estimate: $\tau_{c}$','Interpreter','latex')
        ylabel('MODIS Estimate: $\tau_{c}$','Interpreter','latex')
        title(['Bands $1\& 7$ - RMS: ',num2str(rms_diff_T17)],'Interpreter','latex')
        set(f, 'Position', [0 0 500 500])
        axis square

        % ---------------------------------------------
        % ----------- Plot Bands 2 and 7 --------------
        % ---------------------------------------------

        % MODIS may have used bands 2 and 7 instead of 1 and 7
        % find the minimum and maximum values to create a y=x line

        % min_est = min(est_T27);
        % min_modis = min(modis_T17);
        %
        % max_est = max(est_T27);
        % max_modis = max(modis_T17);
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
        % plot(est_T27,modis_T17,'m.')
        % xlabel('My Estimate: \tau_{c}')
        % ylabel('MODIS Estimate: \tau_{c}')
        % title(['Bands 2&7 - RMS: ',num2str(rms_diff_T27),' \mum'])
        % set(f, 'Position', [0 0 1000 400])



    elseif inputs.bands2plot(nn,1)==1 && inputs.bands2plot(nn,2)==6



        % ---------------------------------------------
        % ----------- Plot Bands 1 and 6 --------------
        % ---------------------------------------------


        modis_T16 = truth_estimate_table.modisT16;
        modis_T16_uncert = modis_T16.*(truth_estimate_table.modisT16_uncert./100);

        est_T16 = truth_estimate_table.estT16;



        square_diffT16 = truth_estimate_table.squareDiffT16; % the absolute difference between my estimate and the modis estimate
        rms_diff_T16 = sqrt(mean(square_diffT16));

        % find the minimum and maximum values to create a y=x line

        min_est = min(est_T16);
        min_modis = min(modis_T16);

        max_est = max(est_T16);
        max_modis = max(modis_T16);

        min_global = min([min_est,min_modis]);

        max_global = min([max_est,max_modis]);

        x = linspace((0.9 * min_global),(1.1*max_global),150);


        f = figure; plot(x,x,'k-','Linewidth',1)
        hold on; grid on; grid minor
        errorbar(est_T16,modis_T16, modis_T16_uncert,'vertical','m.','MarkerSize',10)
        xlabel('My Estimate: \tau_{c}')
        ylabel('MODIS Estimate: \tau_{c}')
        title(['Bands 1&6 - RMS: ',num2str(rms_diff_T16),' \mum'])
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
%     [~, index_2find(ii)] = max(square_diffT17);
%
%     square_diffT17(index_2find(ii)) = 0; % set it to a value that will never be chosen!
%
% end


% what do I want to look at with the estimates that deviate the most from
% the modis values?

% first lets plot all pixels again, and then highlight the 10 that deviate
% the most. This allows us to see how they stack up with the full set

% f = figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% plot(est_T17,modis_T17,'m.')
% xlabel('My Estimate: \tau_{c}')
% ylabel('MODIS Estimate: \tau_{c}')
% %plot(est_T17(index_2find),modis_T17(index_2find),'c.')
% legend('Perfect Fit','all pixels','Location','best')
% %legend('Perfect Fit','all pixels',[num2str(num2find),' furthest from line'],'Location','best')
% %title(['Mean Abs Difference: ',num2str(avg_abs_diff),' \mum'])
% title(['Bands 1&7 - RMS: ',num2str(rms_diff_T17),' \mum'])
% set(f, 'Position', [0 0 1000 400])


% find and remove values of tau that modis deems to be greater than 80.
% This is the upper limit I set in my look up tables.

% index_80 = truth_estimate_table.modisT17>80;
% est_T17_80 = est_T17(~index_80);
% modis_T17_80 = modis_T17(~index_80);
%
% % redefine the absolute difference since we altered some values up above.
% % Also, ignore any values where the tau is greater than 80
%
% square_diffT17 = truth_estimate_table.squareDiffT17(~index_80); % the absolute difference between my estimate and the modis estimate
% num2find = 10;
% index_2find = zeros(1,num2find);
%
% for ii = 1:num2find
%
%     [~, index_2find(ii)] = max(square_diffT17);
%
%     square_diffT17(index_2find(ii)) = 0; % set it to a value that will never be chosen!
%
% end

% figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% plot(est_T17_80,modis_T17_80,'m.')
% xlabel('My estimate - r_{e} (\mum)')
% ylabel('MODIS estimate - r_{e} (\mum)')
% plot(est_T17_80(index_2find),modis_T17_80(index_2find),'c.')
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
% % extract the modis estimate and my calculation estimates
% modis_T17 = truth_estimate_table.modisT17(index24);
% modis_T17_uncert = modis_T17.*(truth_estimate_table.modisT17_uncert(index24)./100);
%
% est_T17 = truth_estimate_table.estT17(index24);
%
%
% square_diffT17 = truth_estimate_table.squareDiffT17(index24); % the absolute difference between my estimate and the modis estimate
% rms_diff_T17 = sqrt(mean(square_diffT17));



% ---------------------------------------------
% ----------- Plot Bands 1 and 7 --------------
% ---------------------------------------------

% find the minimum and maximum values to create a y=x line

% min_est = min(est_T17);
% min_modis = min(modis_T17);
%
% max_est = max(est_T17);
% max_modis = max(modis_T17);
%
% min_global = min([min_est,min_modis]);
%
% max_global = min([max_est,max_modis]);
%
% x = linspace((0.9 * min_global),(1.1*max_global),150);


% f = figure; plot(x,x,'k-','Linewidth',1)
% hold on; grid on; grid minor
% errorbar(est_T17,modis_T17, modis_T17_uncert,'vertical','m.','MarkerSize',10)
% xlabel('My Estimate: \tau_{c}')
% ylabel('MODIS Estimate: \tau_{c}')
% title(['Bands 1&7 - RMS: ',num2str(rms_diff_T17)])
% set(f, 'Position', [0 0 1000 400])





end