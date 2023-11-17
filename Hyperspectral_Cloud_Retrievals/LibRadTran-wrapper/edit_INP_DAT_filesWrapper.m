%% ----- Script Wrapper for editing .INP and .DAT files ----



% By Andrew J. Buggee
%% --- EDIT .INP Files ---

% wavelengthBand ='0_65';

% oldFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/',...
%             'LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/',wavelengthBand,'mu_rev2/'];

oldFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/MODIS_08_25_2021/'];



newFolder = oldFolder;





% introduce change variables in the file name
re = 4:4:20; % effective radius
tau_c = [1,5:5:30]; % cloud optical depth

pixel_row = 1000;
pixel_col = 680;
bands2run = [1:7];

[inpNames2] = getMODIS_INPnames_withClouds(solar,pixel_row,pixel_col,bands2run); % names of files to edit

inpNames2 = reshape(inpNames2,1,[]);


%%

%% ---- Edit and Save new Files ----

% step through each file, edit it and save it as a new file



for kk = 1:length(inpNames2)
    
    
    % redefine the old file each time
    oldFile = inpNames2{kk};
    newFile = oldFile;
    
    oldExpr = 'umu 0.0';
    newExpr = 'umu 1.0';
    
    
    edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr);
    
end

%% ---- Edit and Save new Files ----

% step through each file, edit it and save it as a new file



for kk = 1:length(inpNames2)
    
    for ii = 1:length(re)
        for jj = 1:length(tau_c)
            
            % redefine the old file each time
            oldFile = inpNames2{kk};
            newFile = [inpNames2{kk}(1:end-4),'_r_',num2str(re(ii)),'_T_',num2str(tau_c(jj)),'.INP'];
            
            oldExpr = 'umu 0.0';
            
            if re(ii)<10 && tau_c(jj)<10
                newExpr = ['wc_file 1D ../data/wc/WC_r0',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
                
            elseif re(ii)>=10 && tau_c(jj)<10
                newExpr = ['wc_file 1D ../data/wc/WC_r',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            elseif re(ii)>=10 && tau_c(jj)>=10
                newExpr = ['wc_file 1D ../data/wc/WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            elseif re(ii)<10 && tau_c(jj)>=10
                newExpr = ['wc_file 1D ../data/wc/WC_r0',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
                
            end
            
            
            
            edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr);
            
        end
    end
    
end

%% --------------------------------------------------------
% ----- TAKE A BREAK MAN!!! --------------------



%%

%% ---- Edit and Save new Files ----

% step through each file, edit it and save it as a new file

for ii = 1:length(re)
    for jj = 1:length(tau_c)
        
        % redefine the old file each time
        oldFile = ['WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
        
        oldExpr = ['wc_file 1D ../data/wc/WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj))];
        
        if re(ii)<10 && tau_c(jj)<10
            
            newFile = ['WC_r0',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            newExpr = ['wc_file 1D ../data/wc/WC_r0',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            
        elseif re(ii)>=10 && tau_c(jj)<10
            
            newFile = ['WC_r',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            newExpr = ['wc_file 1D ../data/wc/WC_r',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            
        elseif re(ii)>=10 && tau_c(jj)>=10
            
            newFile = ['WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            newExpr = ['wc_file 1D ../data/wc/WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            
        elseif re(ii)<10 && tau_c(jj)>=10
            
            newFile = ['WC_r0',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            newExpr = ['wc_file 1D ../data/wc/WC_r0',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            
        end
        
        
        
        edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr);
        
    end
end








%% --- EDIT .DAT Files ---


oldFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/',...
    'Hyperspectral-Cloud-Droplet-Retrieval-Research/LibRadTran/libRadtran-2.0.4/data/wc/'];

newFolder = oldFolder;






%% ---- Edit and Save new Files ----

oldFolder = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval-Research/',...
    'LibRadTran/libRadtran-2.0.4/data/wc/'];

newFolder = oldFolder;

% introduce change variables in the file name
re = [30]; % effective radius
tau_c = [50]; % cloud optical depth

den = 10^6; % g/m^3 - density of liquid water
dz = 10^3; % meters - depth of the cloud

% step through each file, edit it and save it as a new file

for ii = 1:length(re)
    for jj = 1:length(tau_c)
        
        % redefine the old file each time
        %         oldFile = ['WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
        oldFile = 'file2Edit';
        
        if re(ii)<10 && tau_c(jj)<10
            newFile = ['WC_r0',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
            
        elseif re(ii)>=10 && tau_c(jj)<10
            newFile = ['WC_r',num2str(re(ii)),'_T0',num2str(tau_c(jj)),'.DAT'];
        elseif re(ii)>=10 && tau_c(jj)>=10
            newFile = ['WC_r',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
        elseif re(ii)<10 && tau_c(jj)>=10
            newFile = ['WC_r0',num2str(re(ii)),'_T',num2str(tau_c(jj)),'.DAT'];
            
        end
        
        % Define old expression to edit
        
        %         oldExpr = tau_c(jj)*2*(re(ii)*10^(-6))*den/(dz); % I forgot the 3 in the old Expr
        oldExpr = '      2.000   1.0000  2'; % I forgot the 3 in the old Expr
        
        % define new value for LWC
        
        newExpr2 = tau_c(jj)*2*(re(ii)*10^(-6))*den/(3*dz);
        
        % we need to make sure the new string length is identical to the
        % old one
        if length(num2str(re(ii)))==1
            num = round(newExpr2,4);
            a = floor(num);
            b = num-a; % this tells us if our number is a fraction, or a whole number
            if b==0
                newExpr2 = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
            else
                
                newExpr2 = num2str(num);
            end
            
            while length(newExpr2)<6
                newExpr2 = [newExpr2,'0'];
            end
            
        elseif length(num2str(re(ii)))==2
            num = round(newExpr2,3);
            a = floor(num);
            b = num-a; % this tells us if our number is a fraction, or a whole number
            if b==0
                newExpr2 = [num2str(num),'.0']; % if true, then there is no deicaml point in our string
            else
                
                newExpr2 = num2str(num);
            end
            
            while length(newExpr2)<5
                newExpr2 = [newExpr2,'0'];
            end
            
        end
        
        % sometimes there are still not enough characters, num2str drops
        % zeros at the end
        
        newExpr = ['      2.000   ',num2str(newExpr2),'  ',num2str(re(ii))];
        
        edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr);
        
    end
end


