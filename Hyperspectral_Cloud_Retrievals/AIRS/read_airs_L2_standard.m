%% View the AIRS water vapor profile retrieval

% By Andrew John Buggee

%% Define the file name and grab file info

airs_data_path = '/Users/andrewbuggee/Documents/MATLAB/Matlab-Research/Hyperspectral_Cloud_Retrievals/AIRS/';


filename = 'AIRS.2024.10.02.234.L2.RetStd_IR.v6.0.34.0.G24277154940.hdf';


% retrieve the info hdf info structure
info = hdfinfo([airs_data_path, filename]);

%%

