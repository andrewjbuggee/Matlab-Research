%% ----- Surface Plots for Model and Observations -----



% By Andrew J. Buggee

%%

function surfPlots4modisModel_andObs(X,Y,modelData,observation_data,band2Plot)


figure; s1 = surf(X,Y,modelData);
xlabel('Optical Depth')
ylabel('Effective Radius (\mum)')
zlabel(['Reflectance (\lambda = ',num2str(band2Plot),' nm)'])
s1.EdgeColor = 'k';
s1.EdgeAlpha = 0.2;
colorbar


% now plot the reflectance value for the observation

hold on
s2 = surf(X,Y,observation_data);
s2.EdgeColor = 'k';
s2.EdgeAlpha = 0;
set(gcf,"Position", [0 0 1000 700])

legend('Model Reflectance','Measured Reflectance','Location','best')




end