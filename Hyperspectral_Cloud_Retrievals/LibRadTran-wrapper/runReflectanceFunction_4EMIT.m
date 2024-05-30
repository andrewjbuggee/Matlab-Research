%% ----- Run the Reflectance Function calculation over a range of re and tau -----


% the input names file must chagne tau across the column space, and must
% change r across row space

% By Andrew J. Buggee

%%

function [R,Rl] = runReflectanceFunction_4EMIT(inputs, names, spectral_response)

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

save_calculated_reflectances_folder = inputs.folder2save.reflectance_calcs;
save_calculated_reflectances_filename = inputs.reflectance_calculations_fileName;

inputFileNames = names.inp;
outputFileNames = names.out;

% We only need the spectral response functions for the bands being using in
% this analysis
spectral_response_2run = cell(1, length(inputs.bands2run));

for nn = 1:length(inputs.bands2run)

    spectral_response_2run{nn} = spectral_response{inputs.bands2run(nn)};

end

length_tau = size(outputFileNames,3);


R = zeros(size(inputFileNames)); % each value here is integrated over the band provided
Rl = cell(size(inputFileNames)); % each value here is the spectral reflectance over the entire band

% Don't compute refelctivity with uvSpec
computeReflectivity = false;


tic
% first step through pixel space
for pp = 1:size(inputFileNames,1)

    % next step through the band dimension
    for bb = 1:size(inputFileNames,4)


        % next step through the different values of effective radius per band
        parfor rr = 1:size(inputFileNames,2)
        %for rr = 1:size(inputFileNames,2)

            % there is a new geometry setting every time we switch
            % pixels. So we need a new input settings for when we
            % switch pixels

            % --- For now, calculate inputSettings every time ---

            % start by running uvspec for a single pixel, a single band, a
            % single effective radius an for all opticl depths
            if iscell(inputFileNames)==true

                [inputSettings] = runUVSPEC(libRadTran_INP_OUT,inputFileNames{pp,rr,:,bb},outputFileNames{pp,rr,:,bb});
                
            elseif ischar(inputFileNames)==true

                [inputSettings] = runUVSPEC(libRadTran_INP_OUT,inputFileNames(pp,rr,:,bb),outputFileNames(pp,rr,:,bb));
            end

            % read in the data structure and calculate reflectance function

            ds = cell(1,length(outputFileNames(pp,rr,:,bb)));

            % --- IMPORTANT ----
            % for the for loop to run, we must define the range of the
            % forloop using a fixed value outside of the parfor loop. The
            % range cannot be defined by calling a sliced variable

            % Finally, step through the different values of optical thickness
            for tt = 1:length_tau

                [pp,rr,tt,bb]

                % Read the UV spec calculation into an array
                [ds{tt},~,~] = readUVSPEC(libRadTran_INP_OUT,outputFileNames{pp,rr,tt,bb},inputSettings(tt+1,:),...
                    computeReflectivity); % headers don't change per iteration

                % ----------- Compute the reflectance function -----------
                [R(pp,rr,tt,bb),~] = reflectanceFunction_4EMIT(inputSettings(tt+1,:), ds{tt},...
                    spectral_response_2run{bb});


            end


        end



    end

end
toc


% Change to the proper MODIS data folder
cd(save_calculated_reflectances_folder);

% Check to see that the file exists
if isfile(save_calculated_reflectances_filename)==true

    % If this file exists, create a new one using 'num' to signify the
    % number of times this has been run
    nn=2;
    % Create a new file name
    save_calculated_reflectances_filename = [save_calculated_reflectances_filename(1:end-4),...
        '_num',num2str(nn), '.mat'];

    while isfile(save_calculated_reflectances_filename)==true

        nn = nn+1;
        % Create a new file name
        save_calculated_reflectances_filename = [save_calculated_reflectances_filename(1:end-4),...
            '_num',num2str(nn), '.mat'];

    end
end

save(save_calculated_reflectances_filename, "inputs", "R"); % save inputSettings to the same folder as the input and output file




end

