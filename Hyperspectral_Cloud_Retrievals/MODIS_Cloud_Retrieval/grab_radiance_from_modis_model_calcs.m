%% ----- Find Radiance from pre-computed uvspec files that were computed for MODIS bi-spectral algorithm -----


function radiance = grab_radiance_from_modis_model_calcs(folderName,names)

inputFileName = names.inp;
outputFileName = names.out;


%% ---- GRAB INPUT SETTINGS ----





%% --- Find the Radiance values ---

% first step through pixel space
for pp = 1:size(inputFileName,1)
    
    
    
    
    
    
    % next step through the band dimension
    for bb = 1:size(inputFileName,4)
        
        
        % ------ LOAD INPUT SETTINGS ------
        
        inputSettings = cell(1,7);
        
        inputSettings{1,1} = 'Solver Type';
        inputSettings{1,2} = 'Cos(zva)';
        inputSettings{1,3} = 'Azimuthal Angle';
        inputSettings{1,4} = 'Solar Zenith Angle';
        inputSettings{1,5} = 'Solar Azimuthal Angle';
        inputSettings{1,6} = 'Sensor Altitude (km)';
        inputSettings{1,7} = 'Source Wavelength (nm) and Irradiance';
        
        textFile = fileread([folderName,inputFileName{pp,1,1,bb}]); % the change r and tau don't effect the input settings
        
        expr1 = '[^\n]*rte_solver [^\n]*';
        expr2 = '[^\n]*umu [^\n]*';
        expr3 = '[^\n]*phi [^\n]*';
        expr4 = '[^\n]*sza [^\n]*';
        expr5 = '[^\n]*phi0 [^\n]*';
        expr6 = '[^\n]*zout [^\n]*';
        expr7 = '[^\n]*source [^\n]*';
        expr8 = '[^\n]*wavelength [^\n]*';
        
        match1 = regexp(textFile,expr1,'match'); % find rte_solver typ
        match2 = regexp(textFile,expr2,'match'); % find consine of viewing angle vector
        match3 = regexp(textFile,expr3,'match'); % find azimuth viewing angle vector
        match4 = regexp(textFile,expr4,'match'); % find the solar zenith angle
        match5 = regexp(textFile,expr5,'match'); % find the solar azimuth angle
        match6 = regexp(textFile,expr6,'match'); % find the sensor altitude
        match7 = regexp(textFile,expr7,'match'); % find the source file
        match8 = regexp(textFile,expr8,'match'); % find the wavelength range in order to trim the source file
        
        index1_space1 = regexp(match1{1},'\s[a-z]+'); % find the spaces
        index1_space2 = regexp(match1{1},'[a-z]\s+'); % the second index should be the last letter in the solver type
        
        index2_space1 = regexp(match2{1},'\s[0123456789]+'); % find a space that is followed by a number
        index2_space2 = regexp(match2{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index3_space1 = regexp(match3{1},'\s[0123456789]+'); % find a space that is followed by a number
        index3_space2 = regexp(match3{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index4_space1 = regexp(match4{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index4_space2 = regexp(match4{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index5_space1 = regexp(match5{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index5_space2 = regexp(match5{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index6_space1 = regexp(match6{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index6_space2 = regexp(match6{1},'[0123456789]\s+'); % find a space that comes after a number
        
        index7_space1 = regexp(match7{1},'\s[a-z]'); % find the spaces
        index7_space2 = regexp(match7{1},'[a-z]\s'); % Brackets treat the symbol literally. number of decimals tells us how many values there are in the vector
        index7_file1 = regexp(match7{1},'flux[/][a-z]'); % find the locaition a letter follows two dots and a forward slash
        index7_file2 = regexp(match7{1},'[.]dat');
        
        index8_space1 = regexp(match8{1},'\s[0123456789]+'); % There is only 1 value for the solar zenith angle
        index8_space2 = regexp(match8{1},'[0123456789]\s+'); % find a space that comes after a number
        
        
        % determine the rte_solver type
        rte_solver = match1{1}(index1_space1(1)+1:index1_space2(2));
        
        % determine the umu vector
        umuStr = cell(1,length(index2_space1));
        
        for ii = 1:length(index2_space1)
            umuStr{ii} = match2{1}(index2_space1(ii)+1:index2_space2(ii));
        end
        
        umuVec = str2double(umuStr);
        
        % determine the phi vector
        phiStr = cell(1,length(index3_space1));
        
        for ii = 1:length(index3_space1)
            phiStr{ii} = match3{1}(index3_space1(ii)+1:index3_space2(ii));
        end
        
        phiVec = str2double(phiStr);
        
        % find the solar zenith angle
        sza = match4{1}(index4_space1+1:index4_space2);
        sza = str2double(sza);
        
        % find the solar azimuth angle
        saz = match5{1}(index5_space1+1:index5_space2(2));
        saz = str2double(saz);
        
        % find the sensor altitude
        
        
        if isempty(index6_space1)==false
            zout = match6{1}(index6_space1(1)+1:index6_space2(2)-1); % this would be for a numeric value
            zout = str2double(zout);
        elseif isempty(index6_space1)==true
            indexString1 = regexp(match6{1},'\s[a-z][a-z][a-z]'); % this would be for a string value like toa or boa
            indexString2 = regexp(match6{1},'[a-z]\s');
            zout_str = match6{1}(indexString1(1)+1:indexString2(2));
            
            if strcmp(zout_str,'toa')==true
                zout = 100;
            elseif strcmp(zout_str,'boa')==true
                zout = 0;
            end
        end
        
        % find the wavelength range of the output file
        wavelength_str= cell(1,length(index8_space1));
        
        for ii = 1:length(index8_space1)
            wavelength_str{ii} = match8{1}(index8_space1(ii)+1:index8_space2(ii));
        end
        wavelength = str2double(wavelength_str);
        
        % determine the solar or thermal source file
        % we want to store the source flux as a vector
        % all solar source files will be located in the folder: /Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/solar_flux
        % all thermal source files will be located in the foler:
        if strcmp('solar',match7{1}(index7_space1(1)+1:index7_space2(2)))
            % find which computer you're running on
            if strcmp(folderName(1:15),'/Users/anbu8374')
                folderSolar = ['/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/'];
                
            elseif strcmp(folderName(1:19),'/Users/andrewbuggee')
                folderSolar = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/',...
                    'LibRadTran/libRadtran-2.0.4/data/solar_flux/'];
            else
                error('Dont reconginze the libradtran solar flux folder')
            end
            fileSolar = match7{1}(index7_file1(1)+5:index7_file2(1)+3);
            sourceFile = fileread([folderSolar,fileSolar]);
            
            exprSource = '[^\n]*[\d][\d] [^\n]*'; % look for the new lines with atleast two digits in a row
            matchSource = regexp(sourceFile,exprSource,'match'); % find rte_solver typ
            
            source = zeros(length(matchSource),2);
            for ii = 1:length(matchSource)
                source(ii,:) = str2num(matchSource{ii});
            end
            
            % now we clip source to match the length of our wavelength vector
            % if we run a monochromatic calculation, we do the following first.
            % Then, for multispectral calculations
            if length(wavelength)==1
                indexSource = source(:,1)==round(wavelength); % can only have integer values for wavelength
                source = source(repmat(indexSource,1,2));
                source = reshape(source,size(source,1)/2,[]);
                
            elseif length(wavelength)>1
                
                indexSource = source(:,1)>=wavelength(1) & source(:,1)<=wavelength(2);
                source = source(repmat(indexSource,1,2));
                source = reshape(source,size(source,1)/2,[]);
            end
        end
        
        % Pull all input settings into a cell array
        inputSettings{1,1} = rte_solver;
        inputSettings{1,2} = umuVec;
        inputSettings{1,3} = phiVec;
        inputSettings{1,4} = sza;
        inputSettings{1,5} = saz;
        inputSettings{1,6} = zout;
        inputSettings{1,7} = source;
        
        % --------------------------------------------------
        
        
        % next step through the different values for r per band
        for rr = 1:size(inputFileName,2)
            
            
            
            
            
            % Finally, step through the different values of optical thickness
            for tt = 1:length_tau
                
                [pp,rr,tt,bb]
                
                [ds,~,~] = readUVSPEC(inp_folder,outputFileName{pp,rr,tt,bb},inputSettings(tt+1,:)); % headers don't change per iteration
                
            end
            
            
        end
        
        
        
    end
    
end




end

