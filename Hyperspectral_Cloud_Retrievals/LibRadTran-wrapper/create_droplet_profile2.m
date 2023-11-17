%% This function will create a cloud droplet profile based on one of four
% physical assumptions. These are taken from S.Platnick's 2000 paper

% INPUTS:
%   (1) re_topBottom - effective droplet radius (microns) - these are the
%   boundary values for the profile you wish to create. You enter them as
%   a vector in the order specified by the variable name [re_top,
%   re_bottom].

%   (2) zT - the vertical independent variable, which is 
%   either geometric altitude or optical depth. If using geometric altitude
%   (z), this variable should be defined in units of kilometers and should
%   start with z_bottom and progress towards z_top. If using optical depth
%   (tau), this vector should start with 0 (cloud top) and end with the
%   total cloud optical thickness.

%   (3) independent varaible - a string that tells the code what the
%   vertical indepdendent variable is. There are two possible inputs here:
%       (a) 'altitude' - tells the code to compute r as a function of
%       geometric height
%       (b) 'optical_depth' - tells the code to compute r as a function of
%       optical thickness within the cloud

%   (4) constraint - the physical constraint (string) - there are four
%   different string options for a physical constraint:
%       (a) 'adiabatic' - this assumption forces the liquid water content to
%       be proportionl to z, the altitude.
%       (b) 'subadiabatic_aloft' - this assumption assumes there is
%       increasing entrainment and drying towards the cloud top.
%       (c) 'linear_with_z' - this constraint forces the effective droplet profile
%       to behave linearly with z (re(z)~z). Physically we are forcing subadiabtatic
%       behavior at mid-levels.
%       (d) 'linear_with_tau' - this constraint forces the effective
%       droplet radius to have linearly with optical depth (re(z)~tau).
%       Physically, this too forces subadiabatic behavior at mid-levels.

% OUTPUTS:
%   (1) re - effective droplet radius profile (microns)


% By Andrew John Buggee
%%

function re = create_droplet_profile2(re_topBottom, zT, independentVariable, constraint)


% ------------------------------------------------------------
% ---------------------- CHECK INPUTS ------------------------
% ------------------------------------------------------------

% Check to make sure there are 4 inputs, droplet radius, cloud optical
% depth, and the altitude vector associated with this cloud


if nargin~=4
    error([newline,'Not enough inputs. Need 4: droplet effective radius at cloud top and bottom,',...
        'the independent variable, a string describing which independent variable is used,',...
        ' and thermodynamic constraint.', newline])
end

% Check to make sure re values are greater than 0

if any(re_topBottom<0)
    
    error([newline,'Effective Radius must be greater than 0.', newline])
end


% Check to make sure the thermodynamic constraint is one 1 four
% possibilites

if strcmp(constraint, 'adiabatic')==false && strcmp(constraint, 'subadiabatic_aloft')==false && strcmp(constraint, 'linear_with_z')==false && strcmp(constraint, 'linear_with_tau')==false
    
    error([newline,'I dont recognize the constraint.', newline])
end


% Check to make sure the independent variable is one of two possibilities

if strcmp(independentVariable, 'altitude')==false && strcmp(independentVariable, 'optical_depth')==false 
    
    error([newline,'I dont recognize the independent variable.', newline])
end




%%



% The boundary values are altered using a powerlaw to get the constraint we
% want

% boundary conditions for r as a function of tau

a0 = @(x) re_topBottom(1)^(2 + 3/x);
a1 = @(x) re_topBottom(1)^(2 + 3/x) - re_topBottom(2)^(2 + 3/x);

% boundary conditions for r as a function of z

b0 = @(x) re_topBottom(2)^(3/x);
b1 = @(x) re_topBottom(1)^(3/x) - re_topBottom(2)^(3/x);

if strcmp(constraint,'subadiabatic_aloft')
    
    % if the profile chosen is subadiabatic aloft, then 0<x<1
    x = 1/2;
    
    
    if strcmp(independentVariable,'optical_depth')==true
        re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));
        
    elseif strcmp(independentVariable,'altitude')==true
        
        re = (b0(x) + b1(x) * (zT - zT(1))./(zT(end) - zT(1))).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
    end
    
    
elseif strcmp(constraint,'adiabatic')
    
    % if the profile chosen is adiabatic, then x=1
    x = 1;
    
    
    if strcmp(independentVariable,'optical_depth')==true
        re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));
        
    elseif strcmp(independentVariable,'altitude')==true
        
        re = (b0(x) + b1(x) * (zT - zT(1))./(zT(end) - zT(1))).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
    end
    
    
elseif strcmp(constraint,'linear_with_z')
    
    % linear with z needs an x of 3
    x = 3;
    
    if strcmp(independentVariable,'optical_depth')==true
        re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));
        
    elseif strcmp(independentVariable,'altitude')==true
        
        re = (b0(x) + b1(x) * (zT - zT(1))./(zT(end) - zT(1))).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
    end
    
    
elseif strcmp(constraint,'linear_with_tau')
    
    % linear with z needs an x of -3
    x = -3;
    
    if strcmp(independentVariable,'optical_depth')==true
        re = (a0(x) - a1(x)*(zT./zT(end))).^(x/(2*x + 3));
        
    elseif strcmp(independentVariable,'altitude')==true
        
        re = (b0(x) + b1(x) * (zT - zT(1))./(zT(end) - zT(1))).^(x/3);                      % droplet profile in geometric coordniate system for an adiabatic cloud
    end
    
    
else
    
    error('I dont recognize the droplet profile you want!')
    
end




end
