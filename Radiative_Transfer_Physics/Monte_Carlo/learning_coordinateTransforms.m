% Coordinate Transformations


clear variables


% Define a set of unit axes
i0 = [1, 0];                    % unit vector along the x axis
j0 = [0, 1];                    % unit vector along the y axis

% Define the function that determines the x and y components of a vector
% due to some angle T with respect to j0 and some vector length L
v = @(L,T) L.*[sind(T), cosd(T)];


%% Test with some vectors

% Define an angle with respect to j0 y axis
T1 = 60;

% Define a vector length
L1 = 1;

% compute the x and y components with respect to the original axes
v1 = v(L1,T1);


% Make a plot of the first vector
figure; plot([0,v1(1)],[0,v1(2)]); grid on; grid minor
xlim([-2,4])
ylim([-2,4])
xline(0,'LineWidth',0.5,'Color','black')
yline(0,'LineWidth',0.5,'Color','black')

%% Compute v_sum using matrix transformations



% Create a new ste of unit vectors so that the new j vector points along v
% The new unit vectors relate to the old unit vectors in the following way
M_cc = @(T) [cosd(T), -sind(T);...
            sind(T), cosd(T)];

M_c = @(T) [cosd(T),  sind(T);...
            -sind(T), cosd(T)];



%% Define the new coordinate system

j1 = M_c(T1)*j0';
i1 = M_c(T1)*i0';


hold on
plot([0, j1(1)], [0, j1(2)],':')
plot([0, i1(1)], [0, i1(2)],':')

%% Compute the second vector using the new coordinate transformations

% Make sure all angles are defined to be clockwise from the y axis
angle_offset_shift = 0;

% Define an angle with respect to the new j1 axis
T2 = 45 + angle_offset_shift;

% Define a vector length
L2= 1;

% compute the x and y components with respect to the new coordinate system
v2 = v(L2,T2);
v2_newCoordinates = M_c(T1)*v(L2,T2)';



% plot the next vector
hold on;
plot([0,v2_newCoordinates(1)],[0, v2_newCoordinates(2)])


% Add the vectors together and plot them

% compute the sum of the two vectors
v_sum = v1 + v2_newCoordinates';

% plot the next vector
hold on;
plot([v1(1),v_sum(1)],[v1(2), v_sum(2)])


%% Try a third vector! Do I need two matrices?


% Make sure all angles are defined to be clockwise from the y axis
angle_offset_shift = 0;

% Define an angle with respect to the new j1 axis
T3 = 45 + angle_offset_shift;

% Define a vector length
L3 = 1;

% compute the x and y components with respect to the new coordinate system
v3 = v(L3,T3);
v3_newCoordinates = M_c(T1)*M_c(T2)*v(L3,T3)';



% plot the next vector
hold on;
plot([0,v3_newCoordinates(1)],[0, v3_newCoordinates(2)])


% Add the vectors together and plot them

% compute the sum of the two vectors
v_sum = [v_sum; v_sum + v3_newCoordinates'];

% plot the next vector
hold on;
plot([v_sum(1,1),v_sum(2,1)],[v_sum(1,2), v_sum(2,2)])



%% Let's try a simulation

% Define a random variable that gives the length of each vector
tau = @(x) -log(1 - x);

% Define a random variable that gives mu=cos(theta) at the end of each
% vector
mu = @(g, x) 1/(2*g) * (1 + g^2 - ((1 - g^2)./(1 - g + 2*g*x)).^2);


% Define the number of vectors you want to plot
N_vectors = 2;

% Plot them!

V_positions = zeros(N_vectors, 2);
tau_draw = zeros(1, N_vectors);
mu_draw = zeros(1, N_vectors);





figure;
xline(0,'LineWidth',0.5,'Color','black')
yline(0,'LineWidth',0.5,'Color','black')


for nn = 1:N_vectors

    %tau_draw(nn) = tau(rand(1));
    tau_draw(nn) = 1;
    %mu_draw(nn) = mu(0.001, rand(1));


    if nn==1

        V_positions(nn,:) = tau_draw(nn).*[sqrt(1 - mu_draw(nn)^2), mu_draw(nn)];

        plot([0, V_positions(nn,1)], [0, V_positions(nn,2)])
        hold on

    else


        V_positions(nn,:) = V_positions(nn-1,:) + tau_draw(nn).*[sqrt(1 - mu_draw(nn)^2), mu_draw(nn)];

        plot([V_positions(nn-1,1), V_positions(nn,1)], [V_positions(nn-1,2), V_positions(nn,2)])
        hold on

    end


end

grid on; grid minor


%% Let's try a simulation using i and j coordniate transformation

clear variables


% Define the function that determines the x and y components of a vector
% due to some angle T with respect to j0 and some vector length L

%v = @(L,Mu) L.*[sqrt(1 - Mu^2), Mu];
v = @(L,Mu) L.*[rand_plusORminus_one(1,1)*sqrt(1 - Mu^2), Mu];



% clockwise rotation
M_c = @(Mu) [Mu, sqrt(1 - Mu^2);...
             -sqrt(1 - Mu^2), Mu];

% counter-clockwise rotation
M_cc = @(Mu) [Mu, -sqrt(1 - Mu^2);...
             sqrt(1 - Mu^2), Mu];


% Define a random variable that gives the length of each vector
tau = @(x) -log(1 - x);

% Define a random variable that gives mu=cos(theta) at the end of each
% vector
mu = @(g, x) 1/(2*g) * (1 + g^2 - ((1 - g^2)./(1 - g + 2*g*x)).^2);


% Define the number of vectors you want to plot
N_vectors = 100;

% Plot them!

M_transformation = eye(2);
position_in_new_coordinates = zeros(N_vectors, 2);
V_positions = zeros(N_vectors, 2);
tau_draw = zeros(1, N_vectors);
mu_draw = zeros(1, N_vectors);
%mu_draw = cos([pi/2, pi/2, pi/2, pi/2]);


figure;


for nn = 1:N_vectors

    %tau_draw(nn) = tau(rand(1));
    tau_draw(nn) = 1;
    mu_draw(nn) = mu(0.85, rand(1));

    % Compute the position in the new coordinate system with respect to the
    % previous vector
    position_in_new_coordinates(nn,:) = v(tau_draw(nn), mu_draw(nn));

    if nn == 1

        % compute the x and y components with respect to the new coordinate system
        V_positions(nn,:) = position_in_new_coordinates(nn,:);
        


        % Plot 
        plot([0, V_positions(nn,1)], [0, V_positions(nn,2)])
        hold on
        drawnow
        pause(0.5)

    else
        
       

%         if nn ==2 || nn==4
%             position_in_new_coordinates(nn,1) = -1*position_in_new_coordinates(nn,1);
%         end



        if position_in_new_coordinates(nn-1,1)<0
            % If x is less than 0, we perform a counter-clockwise transformation
            M_transformation = M_transformation * M_cc(mu_draw(nn-1));
        else
            % If x is greater than 0, we perform a clockwise transformation 
            M_transformation = M_transformation * M_c(mu_draw(nn-1));
        end


        V_positions(nn,:) = V_positions(nn-1,:) + (M_transformation * position_in_new_coordinates(nn,:)')';
        


        % plot
        plot([V_positions(nn-1,1), V_positions(nn,1)], [V_positions(nn-1,2), V_positions(nn,2)])
        hold on
        drawnow
        pause(0.5)


    end



end
xline(0,'LineWidth',0.5,'Color','black')
yline(0,'LineWidth',0.5,'Color','black')
% ylim([-3,3])
% xlim([-3,3])
grid on; grid minor



