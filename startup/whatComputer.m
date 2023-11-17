%% Function to determine which computer is running uvspec


% By Andrew J. Buggee

%%

function [username] = whatComputer()

[status,username] = system('whoami');

username = username(1:end-1);

% determine if we're on the CURC supercomputer or my LASP mac desktop
if strcmp(username, 'anbu8374')==true && strcmp(matlabroot, '/Applications/MATLAB_R2022b.app')==true

    % Then we're using my Mac desktop at LASP
    % Keep this username!

elseif strcmp(username, 'anbu8374')==true && strcmp(matlabroot, '/scratch/local/MATLAB')==true

    % Then were on the super computer! Change the username to reflect this
    username = 'curc';

end



if status ~= 0
    error(['Status of command returned value of ',num2str(status)])
end


end

