function [] = plot_rTop_rBot_iterative_progress(GN_outputs, pixels2use, modis)


for ii = 1:length(pixels2use.res1km.linearIndex)

    f = figure; subplot(1,2,1)

    plot(GN_outputs.rms_residual(:,ii),'.-', 'LineWidth',2);
    grid on; grid minor;
    title(['Pixel ',num2str(ii)], 'Interpreter','latex')
    xlabel('Iteration', 'Interpreter','latex')
    ylabel('RMS Residual', 'Interpreter','latex')

    subplot(1,2,2)
    plot(GN_outputs.retrieval(1,:,ii)','.-', 'LineWidth',2);
    hold on; grid on; grid minor
    plot(GN_outputs.retrieval(2,:,ii)','.-', 'LineWidth',2);
    xlabel('Iteration', 'Interpreter','latex')
    ylabel('$r_e$ ($\mu m$)', 'Interpreter','latex')

    yyaxis right
    plot(GN_outputs.retrieval(3,:,ii)','.-', 'LineWidth',2);
    yline(modis.cloud.optThickness17(pixels2use.res1km.linearIndex(ii)), 'k:','TBLUT $\tau_c$','linewidth',2, 'Interpreter','latex')
    ylabel('$\tau_c$', 'Interpreter','latex')
    yyaxis left
    yline(modis.cloud.effRadius17(pixels2use.res1km.linearIndex(ii)), 'k--','TBLUT $r_e$', 'linewidth',2, 'Interpreter','latex')

    legend('$r_{top}$', '$r_{bot}$','Location','best', 'Interpreter','latex')




    set(f,"Position", [0 0 1000 400])

end



end