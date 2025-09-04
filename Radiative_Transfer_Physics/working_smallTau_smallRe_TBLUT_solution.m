
lgnd_str = cell(1,25);

figure;

for rr = 1:25

    semilogy(ds.wavelength, ds.ssa(:,rr), 'Color', mySavedColors(1+(rr-1), 'fixed'))

    hold on

    lgnd_str{rr} = ['$r_{e}$ = ', num2str(ds.r_eff(rr)), ' $\mu m$'];

end

grid on;

grid minor

legend(lgnd_str, 'Interpreter','latex', 'Location','best', 'FontSize', 20)

title('Single Scattering Albedo of liquid water', 'Interpreter','latex', 'FontSize',30)
ylabel('$\varpi$', 'Interpreter','latex', 'FontSize',35)
xlabel('$\lambda$ (nm)', 'Interpreter','latex', 'FontSize',35)

set(gcf,'Position',[0 0 1350 750])


%%

idx_wl_1 = find(ds.wavelength == 500);  % The first wl used in the TBLUT retrieval
idx_wl_2 = find(ds.wavelength == 2131); % the second wl used in the TBLUT retrieval

% step through different tau values
% step through different effective radius values
tau_c = 1:10;
r_e = 1:7;

R = zeros(numel(tau_c), numel(r_e), 2);
T = zeros(numel(tau_c), numel(r_e), 2);
A = zeros(numel(tau_c), numel(r_e), 2);


for tt = 1:numel(tau_c)

    for rr = 1:numel(r_e)

        idx_re = find(r_e(rr)==ds.r_eff);


        [R(tt,rr,1), T(tt,rr,1), A(tt,rr,1)] = two_stream_RT_boundaryValues(tau_c(tt), ds.ssa(idx_wl_1, idx_re),...
            ds.asymParam(idx_wl_1, idx_re), 0);

        [R(tt,rr,2), T(tt,rr,2), A(tt,rr,2)] = two_stream_RT_boundaryValues(tau_c(tt), ds.ssa(idx_wl_2, idx_re),...
            ds.asymParam(idx_wl_2, idx_re), 0);

    end

end

%% plot Nakajima king style!


% Create a legend string
legend_str = cell(1, length(tau_c) + length(r_e));
% set the first length(tau_c)+1 entries to be empty strings
for ss = 1:length(tau_c)+1
    legend_str{ss} = '';
end


figure;

% Step through each effective radius. Each line represents the reflectance
% at a constant droplet size, with varrying optical thickness.
for rr = 1:length(r_e)

    plot(R(:,rr,1), R(:,rr,2),...
        '.-', 'MarkerSize',50,'LineWidth',1.5);
    hold on

    % store legend string for later
    legend_str{length(tau_c) + rr} = ['$r_e = $', num2str(r_e(rr))];

end

% set up color order for each curve
colororder(mySavedColors(1:length(r_e),'fixed'));



% ------ Plot lines of constant optical depth ------
% Now step through each optical depth. Each line represents the reflectance
% at a constant optical depth, with varrying effective radius.
for tt = 1:length(tau_c)

    x = R(tt,:,1);
    y = R(tt,:,2);

    t = plot(x, y, 'LineStyle','--', 'Color','k');

    % add line label on plot
    if tt==1
        text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",25, "FontWeight","bold")
        hold on
    else
        text(0.995*x(end), 0.8*y(end), num2str(tau_c(tt)),'Interpreter','latex',"FontSize",25, "FontWeight","bold")
        hold on
    end

    % place the dotted line below the lines of constant radius
    uistack(t, 'bottom');

end

% Add text to indicate the black lines are lines of constant optical
% thickness
text(1.05*x(end), 0.85*y(end), '$\tau_c$','Interpreter','latex',"FontSize",25, "FontWeight","bold")



% set up plot stuff
grid on; grid minor
xlabel(['Reflectance 500 $nm$'],Interpreter='latex')
ylabel(['Reflectance 2131 $nm$'],Interpreter='latex')



% set the last string entry to be MODIS value

legend(legend_str, 'Interpreter','latex','Location','best' , 'FontSize', 20, 'FontWeight','bold')
title('Simulated HySICS Reflectance for a vertically homogeneous cloud','Interpreter','latex')
set(gcf,"Position", [0 0 1300 800])



