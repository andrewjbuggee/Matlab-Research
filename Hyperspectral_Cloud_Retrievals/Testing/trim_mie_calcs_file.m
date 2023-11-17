%% Read and Edit the pre-computed mie tables. They are too big to open!!

clear variables


foldername = ['/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/',...
    'LibRadTran/libRadtran-2.0.4/Mie_Calculations/'];

filename = 'Mie_calcs_monodispersed_Wiscombe.OUT';

wavelength_range = 100:3000;                % nanometers
radius_range = 1:500;                       % microns

% After 300 microns, lets jump every 5 microns
desired_r_range = [1:300, 305:5:500];       % microns


% read lines that only begin with the desired wavelength and radius and
% write this into a new file

% Open file and split each line into a new cell array
text = regexp(fileread([foldername,filename]),'\n','split');


% Open new txt file 2 write to
new_filename = 'Mie_calcs_monodispersed_Wiscombe_reduced_file.txt';

% Open and create new file
fileID = fopen([foldername,new_filename], 'w');


for ww = 1:length(wavelength_range)
    for rr = 1:length(desired_r_range)


        if length(num2str(wavelength_range(ww)))==3 && length(num2str(desired_r_range(rr)))==1
            find_expr = ['  ',num2str(wavelength_range(ww)),'.000','     ',num2str(desired_r_range(rr)),'.000000'];

        elseif length(num2str(wavelength_range(ww)))==3 && length(num2str(desired_r_range(rr)))==2
            find_expr = ['  ',num2str(wavelength_range(ww)),'.000','    ',num2str(desired_r_range(rr)),'.000000'];

        elseif length(num2str(wavelength_range(ww)))==3 && length(num2str(desired_r_range(rr)))==3
            find_expr = ['  ',num2str(wavelength_range(ww)),'.000','   ',num2str(desired_r_range(rr)),'.000000'];

        elseif length(num2str(wavelength_range(ww)))==4 && length(num2str(desired_r_range(rr)))==1
            find_expr = [' ',num2str(wavelength_range(ww)),'.000','    ',num2str(desired_r_range(rr)),'.000000'];

        elseif length(num2str(wavelength_range(ww)))==4 && length(num2str(desired_r_range(rr)))==2
            find_expr = [' ',num2str(wavelength_range(ww)),'.000','   ',num2str(desired_r_range(rr)),'.000000'];

        elseif length(num2str(wavelength_range(ww)))==4 && length(num2str(desired_r_range(rr)))==3
            find_expr = ['  ',num2str(wavelength_range(ww)),'.000','  ',num2str(desired_r_range(rr)),'.000000'];

        end



        % find location of expression to change
        write_line_2_new_file = find(contains(text,find_expr)); % find the old expression location




        % fprintf writes lines in our text file from top to botom

        % Write in the data
        fprintf(fileID,'%s \n', text{write_line_2_new_file});



    end
end


% close the new file
fclose(fileID);







