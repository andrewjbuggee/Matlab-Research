function data = read_libRadtran_atm_dat_profiles_ver2(filename)
%READ_AFGLUS Reads the AFGLUS .DAT file into a MATLAB matrix, ignoring comments.
%
%   data = READ_AFGLUS(filename) reads the specified .DAT file and returns
%   a numeric matrix. Lines starting with '#' are ignored.

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read all lines into a cell array
    allLines = {};
    tline = fgetl(fid);
    while ischar(tline)
        allLines{end+1, 1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    
    % Filter out empty lines and comment lines
    validLines = {};
    for i = 1:length(allLines)
        trimmedLine = strtrim(allLines{i});
        % Keep line if it's not empty and doesn't start with '#'
        if ~isempty(trimmedLine) && trimmedLine(1) ~= '#'
            validLines{end+1, 1} = trimmedLine;
        end
    end
    
    % Convert to numeric matrix
    numCols = 9;
    data = zeros(length(validLines), numCols);
    
    for i = 1:length(validLines)
        data(i, :) = sscanf(validLines{i}, '%f', numCols);
    end

end