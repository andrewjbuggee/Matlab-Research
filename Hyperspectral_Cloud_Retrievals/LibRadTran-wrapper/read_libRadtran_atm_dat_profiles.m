function data = read_libRadtran_atm_dat_profiles(filename)
%READ_AFGLUS Reads the AFGLUS .DAT file into a MATLAB matrix, ignoring comments.
%
%   data = READ_AFGLUS(filename) reads the specified .DAT file and returns
%   a numeric matrix. Lines starting with '#' are ignored.

    % Read all lines
    allLines = readlines(filename);

    % Remove empty lines and comment lines
    validLines = allLines(~startsWith(strtrim(allLines), "#") & strlength(strtrim(allLines)) > 0);

    % Convert to numeric matrix
    numCols = 9;
    data = zeros(length(validLines), numCols);

    for i = 1:length(validLines)
        data(i, :) = sscanf(validLines(i), '%f', numCols);
    end

end
