%% Simple exponential analytical weighting function
clear variables

% The depth of the cloud is defined between y=0 and y=H
H = 1;  % km

% the weighting function is w(y) = exp(a*H) - 1
% The integral of w(y) = 1
% solve for a
syms a
eqn = exp(a*H)/a - H - 1/a ==1;

% solve for a
solve_a = solve(eqn, a, 'Real', true);

% take the real solution to be the value for a
a = eval(solve_a);

% Test to make sure w(y) integrates to 1
y = linspace(0, H, 1000);           % km
w_sat = exp(a*y)-1;

disp(['The integral of w(y) between o and H is: ', num2str(trapz(y, w_sat) )])

% plot the result
figure; plot(w_sat, y); grid on; grid minor
xlabel('w(y)'); ylabel('y')


%% Solve for an analytical gaussian weighting function for the helicopter

% helicopter moves only in the y direction at a velocity of 5 m/s
vy = 5;     % m/s
T = (H*1000)/vy;        % sec - total time to clear the cloud
t = linspace(0, T, 500);     % seconds it takes for the helicopter to clear the cloud

mu = (vy/1e3)*t;    % km - position of the helicopter as it ascends through a cloud
sig = H/100;    % km - the standard deviation, covering 66% of the weighting function - H/100


% Solve just one verison of the weighting function for the helicopter and
% plot this example
mu_idx = 100;

c = 2 / (sig * sqrt(pi) * (erf((H - mu(mu_idx))/sig) + erf(mu(mu_idx)/sig)));

w_heli = c* exp(-(y - mu(mu_idx)).^2 ./sig^2);

disp(['The integral of w(y) between o and H is: ', num2str(trapz(y, w_heli) )])


% plot the result
figure; plot(w_heli, y); grid on; grid minor
xlabel('w(y)'); ylabel('y')


%% create a simple 1D cloud model with some random noise that changes with time

% effective droplet radius is linear with cloud height
r_bot = 5;      % microns
r_top = 10;         % microns

% there is a random gaussian component that is only a function of time
% the variance of the spread increases linearly with time
sig_0 = 0.1;        % microns
sig_t = 1;          % microns


sig_time = sig_0 + sig_t*(t/T);


%% Let's plot a movie of the droplet radius profile over time

% compute the profile without the random component
re = r_bot + (r_top-r_bot)*(y/H);           % microns

figure; 
for tt = 1:length(t)
    


    % compute the random gaussian component
    % randn is a guassian random variable with a mean of 0 and a varaince
    % of 1
    G = sig_time(tt)*randn(1, length(re));          % microns

    % add the random noise to the re profile above. We keep the same
    % profile shape for all of time but the variance increases with time
    re_tt = re + G;     % microns



    plot(re_tt, y)
    grid on; grid minor
    title(['t = ', num2str(t(tt)), ' sec'])
    xlim([0 20])
    ylabel('y (km)')
    xlabel('r_e(y, t)  (\mum)')
    drawnow
    pause(0.1)


end


%% The helicopter flies vertically, sampling the droplets every second

% Are the units correct?
% I've dont this before in optical depth space. I don't recall having to
% divide by anything. But that time I was using a weighting function in
% optical depth space. This time its in geometrical space

sampling_time = 3;      % sec - duration of sampling time

instrument_time = t(t<=sampling_time);  % first sampling time seconds

instrument_height = (vy*instrument_time')./1e3;         % km

Y = repmat(y, length(instrument_time), 1);

w_heli = exp(-(repmat(y, length(instrument_time), 1) - instrument_height).^2 ./sig^2);     % 1/km
    


% compute the random gaussian component
% randn is a guassian random variable with a mean of 0 and a varaince
% of 1
G = repmat(sig_time(t<=sampling_time)', 1, length(re)) .* randn(length(instrument_time), length(re));          % microns

% add the random noise to the re profile above. We keep the same
% profile shape for all of time but the variance increases with time
re_tt = repmat(re, length(instrument_time), 1) + G;     % microns


% first integrate over y
int_over_y = trapz(y, re_tt .* w_heli,2)

% integrate over time and over y
trapz(instrument_time, int_over_y, 1)

