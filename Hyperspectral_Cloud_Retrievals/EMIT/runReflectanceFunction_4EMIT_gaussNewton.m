%% ----- Run the Reflectance Function calculation over a range of re and tau -----


% the input names file must chagne tau across the column space, and must
% change r across row space

% By Andrew J. Buggee

%%

function [R,Rl] = runReflectanceFunction_4EMIT_gaussNewton(names, inputs, spectral_response)

% what computer are we using?

userName = whatComputer;

if strcmp(userName,'anbu8374')

    libRadTran_path = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4'];

elseif strcmp(userName,'andrewbuggee')

    libRadTran_path = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
        'LibRadTran/libRadtran-2.0.4'];

else
    error('I dont recognize this computer user name')
end


addpath(libRadTran_path);

% ----- extract inputs -----


libRadTran_INP_OUT = inputs.folder2save.libRadTran_INP_OUT; % where the newly created .inp files will be saved



inputFileNames = names.inp;
outputFileNames = names.out;


R = zeros(size(inputFileNames)); % each value here is integrated over the band provided
Rl = cell(size(inputFileNames)); % each value here is the spectral reflectance over the entire band

% Don't compute refelctivity with uvSpec
computeReflectivity = false;


% --- step through the band dimension ---
%parfor bb = 1:size(inputFileNames, 2)
for bb = 1:size(inputFileNames,2)
    % --- For now, calculate inputSettings every time ---

    % run uvSpec across all wavelengths
    [inputSettings] = runUVSPEC(libRadTran_INP_OUT, inputFileNames{bb}, outputFileNames{bb});


    % read in the data structure and calculate reflectance function
    % Read the UV spec calculation into an array
    [ds,~,~] = readUVSPEC(libRadTran_INP_OUT, outputFileNames{bb},inputSettings(2,:),...
        computeReflectivity); % headers don't change per iteration

    % ----------- Compute the reflectance function -----------
    [R(bb),~] = reflectanceFunction_4EMIT(inputSettings(2,:), ds,...
        spectral_response(bb,:)');

end




end

