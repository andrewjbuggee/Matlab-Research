%% ----- Run the Reflectance Function calculation over a range of re and tau -----


% the input names file must chagne tau across the column space, and must
% change r across row space

% By Andrew J. Buggee

%%

function [R,Rl] = runReflectanceFunction(inputs,names, spectral_response)

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


inp_folder = [libRadTran_path,'/',inputs.INP_folderName]; % where the newly created .inp files will be saved
saveCalcs_filename = inputs.saveCalculations_fileName;
saveCalcs_folder = inputs.savedCalculations_folderName;
inputFileNames = names.inp;
outputFileNames = names.out;

length_tau = size(outputFileNames,3);


R = zeros(size(inputFileNames)); % each value here is integrated over the band provided
Rl = cell(size(inputFileNames)); % each value here is the spectral reflectance over the entire band

% Don't compute refelctivity with uvSpec
computeReflectivity = false;



% first step through pixel space
for pp = 1:size(inputFileNames,1)

    % next step through the band dimension
    for bb = 1:size(inputFileNames,4)


        % next step through the different values for r per band
        parfor rr = 1:size(inputFileNames,2)
        %for rr = 1:size(inputFileNames,2)

            % if pp==1 && bb ==1 && rr == 1

            % there is a new geometry setting every time we switch
            % pixels. So we need a new input settings for when we
            % switch pixels

            % --- For now, calculate inputSettings every time ---

            % start by running uvspec
            [inputSettings] = runUVSPEC(inp_folder,inputFileNames(pp,rr,:,bb),outputFileNames(pp,rr,:,bb));

            % save the input settings file. In this case each band will have the
            % same geometry so we only need to save a single settings file per band
            %                 save([inp_folder,'runUVSPEC_settings.mat'],"inputSettings"); % save inputSettings to the same folder as the input and output file

            %             else
            %                 runUVSPEC(inp_folder,inputFileNames(pp,rr,:,bb),outputFileNames(pp,rr,:,bb));
            %             end



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
                [ds{tt},~,~] = readUVSPEC(inp_folder,outputFileNames{pp,rr,tt,bb},inputSettings(tt+1,:), computeReflectivity); % headers don't change per iteration
                
                % Compute the reflectance function
                %[R(pp,rr,tt,bb),Rl{pp,rr,tt,bb}] = reflectanceFunction(inputSettings(tt+1,:),ds{tt}, spectral_response{bb}(:,2));
                [R(pp,rr,tt,bb),Rl{pp,rr,tt,bb}] = reflectanceFunction_4modis(inputSettings(tt+1,:),ds{tt}, spectral_response{bb}(:,2));

            end


        end



    end

end



% Change to the proper MODIS data folder
cd(saveCalcs_folder);

% Check to see that the file exists
if isfile(saveCalcs_filename)==true
    % save the relfectance calculation
    save(saveCalcs_filename,"R",'-append'); % save inputSettings to the same folder as the input and output file
else
    save(saveCalcs_filename, "inputs", "R");
end



end

