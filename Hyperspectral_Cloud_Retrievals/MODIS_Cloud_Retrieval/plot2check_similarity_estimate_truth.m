
function [] = plot2check_similarity_estimate_truth(estimate,truth)





% find the minimum and maximum values to create a y=x line

min_est = min(estimate);
min_modis = min(truth);

max_est = max(estimate);
max_modis = max(estimate);

min_global = min([min_est,min_modis]);

max_global = min([max_est,max_modis]);

x = linspace((0.9 * min_global),(1.1*max_global),150);


figure; plot(x,x,'w-','Linewidth',1)
hold on; grid on; grid minor
plot(estimate,truth,'m.')
xlabel('My estimate - r_{e} (\mum)')
ylabel('MODIS estimate - r_{e} (\mum)') 

end