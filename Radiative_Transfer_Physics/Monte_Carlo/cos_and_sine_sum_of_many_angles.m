%% Compute the cosine of the sum of a set of angles

% There are 3 options for the input units
%   (1) 'radians'
%   (2) 'degrees'
%   (3) 'mu' where mu = cos(theta) - this is a special case used for monte
%   carlo methods in radiative transfer

function [cos_sum, sin_sum] = cos_and_sine_sum_of_many_angles(angles, units)





% Need the length of the vector
N = length(angles);


% Make zero vectors
cos_sum = zeros(1,N);
sin_sum = zeros(1,N);

% Determine which vectors are negative. These angles need to use the
% subtraction formula.
%angle_sign = sign(angles);


if strcmp(units,'radians')==true

    % angles provided are in units of radians


    for nn=1:N

        if nn==1

            cos_sum(nn) = cos(angles(nn));
            sin_sum(nn) = sin(angles(nn));

        elseif nn==2

            cos_sum(nn) = cos_sum_of_2_angles(angles(1:2),'radians');
            sin_sum(nn) = sin_sum_of_2_angles(angles(1:2),'radians');

        else
            cos_sum(nn) = cos_sum(nn-1)*cos(angles(nn)) - sin_sum(nn-1)*sin(angles(nn));
            sin_sum(nn) = sin_sum(nn-1)*cos(angles(nn)) + cos_sum(nn-1)*sin(angles(nn));
        end

    end



elseif strcmp(units,'degrees')==true

    % angles provided are in units of degrees


    for nn=1:N

        if nn==1

            cos_sum(nn) = cosd(angles(nn));
            sin_sum(nn) = sind(angles(nn));

        elseif nn==2

            cos_sum(nn) = cos_sum_of_2_angles(angles(1:2),'degrees');
            sin_sum(nn) = sin_sum_of_2_angles(angles(1:2),'degrees');

        else
            cos_sum(nn) = cos_sum(nn-1)*cosd(angles(nn)) - sin_sum(nn-1)*sind(angles(nn));
            sin_sum(nn) = sin_sum(nn-1)*cosd(angles(nn)) + cos_sum(nn-1)*sind(angles(nn));
        end

    end




elseif strcmp(units,'mu')==true


    % angles provided are as mu=cos(theta)
    mu = angles;


    % Compute nu = sin(theta)

    nu = sqrt(1 - mu.^2);


    for nn=1:N

        if nn==1

            cos_sum(nn) = mu(nn);
            sin_sum(nn) = nu(nn);

        elseif nn==2

            cos_sum(nn) = mu(nn-1)*mu(nn) - nu(nn-1)*nu(nn);
            sin_sum(nn) = nu(nn-1)*mu(nn) + mu(nn-1)*nu(nn);

        else
            cos_sum(nn) = cos_sum(nn-1)*mu(nn) - sin_sum(nn-1)*nu(nn);
            sin_sum(nn) = sin_sum(nn-1)*mu(nn) + cos_sum(nn-1)*nu(nn);
        end

    end



    

else


    error([newline,'I dont recognize the units. The units must be either "radians", "degrees" or',...
        ' "mu". Please specifiy.',newline])

end





end