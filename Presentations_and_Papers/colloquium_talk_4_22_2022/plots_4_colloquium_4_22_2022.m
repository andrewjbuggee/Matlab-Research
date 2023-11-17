%% Plots made for colloqium talk on 4/22/2022


% By Andrew J. Buggee

%% Plot of different droplet profiles according to Platnick
clear variables
scriptPlotting_wht

% Recreating the Platnick Plot on different droplet profiles (2000, fig 1

z = 0:1:100;
idv = 'altitude';
re_tb = [12,5];

re_ad = create_droplet_profile2(re_tb, z, idv, 'adiabatic');
re_z = create_droplet_profile2(re_tb, z, idv, 'linear_with_z');
re_t = create_droplet_profile2(re_tb, z, idv, 'linear_with_tau');

figure; plot(re_ad,z./z(end)); hold on
plot(re_z, z./z(end)); plot(re_t, z./z(end));
%plot(re_t, z./z(end))
xlabel('Effective Radius (\mu m)'); ylabel(['Relative Geometric',newline,' Cloud Thickness (z/H)'])
grid on; grid minor;
%legend('adiabatic','sub-adiabatic','linear w/ z', 'linear w/ \tau','Location','best')
legend('Adiabatic','$r_e \propto z$','$r_e \propto \tau$','Location','best','interpreter','latex',...
    'Fontsize',20)
xlabel('Effective Radius $$(\mu m)$$','Interpreter','latex');
ylabel(['Relative Geometric',newline,' Cloud Thickness $$(z/H)$$'], 'Interpreter','latex')
title('Theoretical Cloud Droplet Profiles','interpreter','latex');

%% Plot of Single Scattering Albedo

clear variables

scriptPlotting_wht

folder_path = '/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/Mie_Calculations/';

filename = 'Mie_Properties_4_AVIRIS_1nm_sampling_monodispersed.OUT';
format_spec = '%f %f %f %f %f %f %f %f';        % 8 columns of data

file_id = fopen([folder_path,filename]);

data_table = textscan(file_id, format_spec);

index_r = [1,10:10:100];

ssa = reshape(data_table{6},100,[]);


wl = reshape(data_table{1},100,[]);


figure;
plot(wl(1,:), ssa(index_r,:))
grid on; grid minor
xlabel('Wavelength (nm)','Interpreter','latex')
ylabel('Signle Scattering Albedo','Interpreter','latex')
legend(strcat(string(index_r),' \mum'),'location','best')
xlim([250, 2500])


%% Creating 2D Gaussian PDF's to show how Bayesian Inverse Theory Works

% -- MODEL PRIOR ---

clear variables

 % how many model parameters are we interested in?
 
 num_model_parameters = 2;
 
 % lets define the mean for each of these
 % if we pretend the two parameters are temperature and pressure...
 mean_m = [8, 25]; % 290 K and 850 mb are the average values we expect
 
 % what do we expect the variance to be in our estimates?
 var_m = [10, 25]; % variance of temp is 10 K and variance of pressure is 5 mb
 
 % we need to define the covariance matrix, which describes how correlated
 % our model parameters are. For now, lets assume they are independent from
 % one another
 
 S_m = diag(var_m); % diagonal covariance matrix, where the off diagonal elements are 0
 
 % Define the number of sample points along each variable
 
 N = 150;
 
 % lets create the space of values that are used to determine the
 % probability
 
 
r_min_max = [0, 25];
T_min_max = [0, 50];

 x = [linspace(r_min_max(1), r_min_max(2),N); linspace(T_min_max(1), T_min_max(2),N)];
 
 [X1,X2] = meshgrid(x(1,:),x(2,:));
 X = [X1(:), X2(:)]; % turn this into a column vector

 P_prior2measurement = mvnpdf(X,mean_m,S_m);
 P_prior2measurement = reshape(P_prior2measurement,length(x(2,:)),length(x(1,:)));
 
 
 figure; s = surf(x(1,:),x(2,:),P_prior2measurement); colorbar
 xlabel('Effective Radius $$(\mu m)$$','Interpreter','latex'); ylabel('Optical Depth','Interpreter','latex');
 zlabel('Probability','Interpreter','latex');
 s.EdgeAlpha = 0.25;
 title('Prior $$P(\vec{m})$$','Interpreter','latex')
 
 % Set axes limits
 xlim([0,20]);
 ylim([0,40]);
 zlim([0,0.08]);
 
 % -----------------------
 % --- MODEL POSTERIOR ---
 % -----------------------
  
 clear variables

 % how many model parameters are we interested in?
 
 num_model_parameters = 2;
 
 % lets define the mean for each of these
 % if we pretend the two parameters are temperature and pressure...
 mean_m = [8, 25]; % 290 K and 850 mb are the average values we expect
 
 % what do we expect the variance to be in our estimates?
 var_m = [1.5, 3]; % variance of temp is 10 K and variance of pressure is 5 mb
 
 % we need to define the covariance matrix, which describes how correlated
 % our model parameters are. For now, lets assume they are independent from
 % one another
 
 S_m = diag(var_m); % diagonal covariance matrix, where the off diagonal elements are 0
 
 % Define the number of sample points along each variable
 
 N = 150;
 
 % lets create the space of values that are used to determine the
 % probability
 
 
r_min_max = [0, 25];
T_min_max = [0, 50];

 x = [linspace(r_min_max(1), r_min_max(2),N); linspace(T_min_max(1), T_min_max(2),N)];
 
 [X1,X2] = meshgrid(x(1,:),x(2,:));
 X = [X1(:), X2(:)]; % turn this into a column vector

 P_prior2measurement = mvnpdf(X,mean_m,S_m);
 P_prior2measurement = reshape(P_prior2measurement,length(x(2,:)),length(x(1,:)));
 
 
 figure; s = surf(x(1,:),x(2,:),P_prior2measurement); colorbar
 xlabel('Effective Radius $$(\mu m)$$','Interpreter','latex'); ylabel('Optical Depth','Interpreter','latex')
 zlabel('Probability','Interpreter','latex');
 s.EdgeAlpha = 0.25;
 title('Posterior $$P(\vec{m}|\vec{d})$$','Interpreter','latex')
     
 
  % Set axes limits
 xlim([0,20]);
 ylim([0,40]);
 zlim([0,0.08]);
     

