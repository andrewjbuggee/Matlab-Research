% Function solving for sine(a+b)

% sin(a+b) = sin(a)cos(b) + cos(a)sin(b)

% specify whether the angles are in degrees or radians



function output = sin_sum_of_2_angles(angles, units)

% ---------------------------------------------------
% --- Check to make sure angles has a length of 2 ---
% ---------------------------------------------------

if length(angles)~=2
    error([newline, '2 angles are required for the angle input!',newline])
end

% ---------------------------------------------------






if strcmp(units,'radians')==true
    
    % angles provided are in units of radians    
    output = sin(angles(1))*cos(angles(2)) + cos(angles(1))*sin(angles(2));



elseif strcmp(units,'degrees')==true

    % angles provided are in units of degrees    


    output = sind(angles(1))*cosd(angles(2)) + cosd(angles(1))*sind(angles(2));


else

    error([newline,'I dont recognize the units. The units must be either "radians" or "degrees".',...
        'Please specifiy.',newline])

end



end