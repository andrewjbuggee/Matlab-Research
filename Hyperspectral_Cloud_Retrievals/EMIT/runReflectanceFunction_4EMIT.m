%% ----- Run the Reflectance Function calculation over a range of re and tau -----


% the input names file must chagne tau across the column space, and must
% change r across row space

% By Andrew J. Buggee

%%

function [R,Rl, inputs] = runReflectanceFunction_4EMIT(inputs, names, spectral_response)



% ----- extract inputs -----


libRadTran_INP_OUT = inputs.folder2save.libRadTran_INP_OUT; % where the newly created .inp files will be saved

save_calculated_reflectances_folder = inputs.folder2save.reflectance_calcs;
save_calculated_reflectances_filename = inputs.reflectance_calculations_fileName;

inputFileNames = names.inp;
outputFileNames = names.out;


length_tau = size(outputFileNames,3);

% Define the spectral response functions for the desired wavelengths
spectral_response_2run = spectral_response(inputs.bands2run, :);


R = zeros(size(inputFileNames)); % each value here is integrated over the band provided
Rl = cell(size(inputFileNames)); % each value here is the spectral reflectance over the entire band

% Don't compute refelctivity with uvSpec
computeReflectivity = false;


tic
% first step through pixel space
for pp = 1:size(inputFileNames,1)

    % next step through the band dimension
    for bb = 1:size(inputFileNames,4)


        % next step through the different values of effective radius
        parfor rr = 1:size(inputFileNames,2)
        %for rr = 1:size(inputFileNames,2)

            % there is a new geometry setting every time we switch
            % pixels. So we need a new input settings for when we
            % switch pixels

            % --- For now, calculate inputSettings every time ---

            % start by running uvspec for a single pixel, a single band, a
            % single effective radius an for all optical depths
            [inputSettings] = runUVSPEC(libRadTran_INP_OUT, inputFileNames(pp,rr,:,bb), outputFileNames(pp,rr,:,bb));



            % read in the data structure and calculate reflectance function

            ds = cell(1,length(outputFileNames(pp,rr,:,bb)));

            % --- IMPORTANT ----
            % for the for loop to run, we must define the range of the
            % forloop using a fixed value outside of the parfor loop. The
            % range cannot be defined by calling a sliced variable

            % Finally, step through the different values of optical thickness
            for tt = 1:length_tau

                disp(['[pp,rr,tt,bb] = [', num2str(pp), ',' num2str(rr), ',', num2str(tt), ',', num2str(bb), ']'])

                % Read the UV spec calculation into an array
                [ds{tt},~,~] = readUVSPEC(libRadTran_INP_OUT,outputFileNames{pp,rr,tt,bb},inputSettings(tt+1,:),...
                    computeReflectivity); % headers don't change per iteration

                % ----------- Compute the reflectance function -----------
                [R(pp,rr,tt,bb),~] = reflectanceFunction_4EMIT(inputSettings(tt+1,:), ds{tt},...
                    spectral_response_2run(bb,:)');


            end


        end



    end

end
toc


% Change to the proper EMIT data folder
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
        save_calculated_reflectances_filename = [save_calculated_reflectances_filename(1:end-5),...
            num2str(nn), '.mat'];

    end
end

% Save the new filename
inputs.reflectance_calculations_fileName = save_calculated_reflectances_filename;

% save the calculated reflectances and the inputs
save(save_calculated_reflectances_filename, "inputs", "R"); % save inputSettings to the same folder as the input and output file




end

