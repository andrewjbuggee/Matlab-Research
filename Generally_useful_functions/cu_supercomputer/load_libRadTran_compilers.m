%% load the proper modules for libRadTran to run on the CU super computer


% By Andrew John Buggee
%%

function load_libRadTran_compilers

%% Load compilers and libararies


% The CU Cupercomputer OnDemand Core desktop uses the RedHat distribution
% of Linux

% If you'd like to know more about a package, try 'module spider <package
% name>

% a successful command will return a status of 0
% an unsuccessful command will return a status of 1


% use the function 'system' to run commands in the terminal window
% change to the base directory


% %cmnd = ['cd ~'];
% cmnd_folder = ['cd /projects/anbu8374/libRadtran-2.0.5'];
% 
% [status] = system(cmnd_folder);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end



% load purge
cmnd = ['ml purge', newline, 'ml gcc/11.2.0', newline, 'ml gsl/2.7', newline, 'ml netcdf/4.8.1'];
[status] = system([cmnd]);
if status ~= 0
    error(['Status returned value of ',num2str(status)])
end

% % load the GNU compiler collection
% % This incluudes gcc, g++, and gfortran
% cmnd = ['ml gcc/11.2.0'];
% [status] = system([cmnd_folder, ' ; ', cmnd]);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% % Load the GNU scientific library
% cmnd = ['ml gsl/2.7'];
% [status] = system([cmnd_folder, ' ; ', cmnd]);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% 
% % Load the netCDF tools and library package
% cmnd = ['ml netcdf/4.8.1'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% 
% 
% % Set path variable
% cmnd = ['export PATH=/projects/$USER/software/libRadtran-2.0.5/bin:$PATH'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end




% % Load the perl programming language
% % This is needed to run the automatic test
% cmnd = ['module load perl/5.36.0'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end



% % Load texlive, a package needed to load postsctipt files for the
% % documentation
% cmnd = ['module load texlive/2021'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end


% % load the intel package
% cmnd = ['module load intel/2022.1.2'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% % load the Intel MPI compiler
% cmnd = ['module load impi/2021.5.0'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end







%% Switch to the libRadTran folder and compile the software

% % switch to the libRadTran folder!1
% cmnd = ['cd libRadtran-2.0.5'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% % Configure libRadTran!!
% % point to the GSL library
% % /curc/sw/install/gsl/2.7/gcc/11.2.0
% % -L/curc/sw/install/gsl/2.7/gcc/11.2.0/lib -lgsl -lgslcblas -lm
% cmnd = ['./configure'];
% [status] = system(cmnd);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end
% 
% 
% 
% 
% % a successful command will return a status of 0
% % an unsuccessful command will return a status of 1
% 
% 
% 
% [status] = system([cmnd1, ' ; ', cmnd2]);
% if status ~= 0
%     error(['Status returned value of ',num2str(status)])
% end



end