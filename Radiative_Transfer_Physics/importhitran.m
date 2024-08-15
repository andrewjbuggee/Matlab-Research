function [ lines ] = importhitran( par_file )
%IMPORTHITRAN Imports the HITRAN database.
%   lines = IMPORTHITRAN(par_file) imports HITRAN from the file par_file.
%
%   Reads the High-resolution Transmission (HITRAN) database from a "par"
%   file. The par file is a large file (~500 MB) containing line parameters
%   such as transmission wavenumber and line intensity. In the 2008 version
%   of the HITRAN database this file is called 'HITRAN08.par'.
%
%   Tip: once the database has been read, save it as a MATLAB MAT file by
%   typing:
%      save('HITRAN08.mat', 'lines')
%
%   To load the MAT file next time:
%      load('HITRAN08.mat')
%
%   It will typically take a few seconds to load the MAT file, compared to
%   minutes or hours to import the HITRAN "par" file.
%
%   Part of this function is based on a MATLAB script called READ_HITRAN2
%   by H. Motteler.
%
%   References:
%      "HITRAN Database Format 2004-present"
%      Harvard-Smithsonian Center for Astrophysics
%      http://www.cfa.harvard.edu/hitran/formats.html
%
%      The HITRAN 2008 molecular spectroscopic database
%      L.S. Rothman and others
%      http://www.cfa.harvard.edu/hitran/docs.html
%
%   Iain Robinson, School of Geosciences, University of Edinburgh
%   4th December 2014
	
% Check that the file (e.g. HITRAN08.par) exists.
if ~exist(par_file, 'file')
    error('The HITRAN database file\n\t%s\nwas not found. This file is called something like HITRAN08.par', par_file);
end

% Open the file.
[fileID, message] = fopen(par_file, 'r', 'native', 'US-ASCII');
if fileID == -1
    close(h)
    error('The HITRAN database file\n\t%s\nexists, but could not be read. The operating system gave the message\n\t%s', par_file, message);
end

% Define some properties of the file.
recordLength = 160;
recordLengthIncludingLineEndings = recordLength + 2; % Lines end with CR+LF.

% Check the length of the file.
fseek(fileID, 0, 'eof');
fileSizeInBytes = ftell(fileID);
fseek(fileID, 0, 'bof');

% Work out the number of lines in the file.
numberOfLines = fileSizeInBytes / recordLengthIncludingLineEndings;

% Preallocate arrays to store the (known) number of lines. As this may take
% some time, show a waitbar on the screen.
h = waitbar(0, sprintf('%d lines in file. Allocating memory...', numberOfLines));
lines.moleculeNumber = zeros(numberOfLines,1);
lines.isotopologueNumber = zeros(numberOfLines,1);
lines.transitionWavenumber = zeros(numberOfLines,1);
lines.lineIntensity = zeros(numberOfLines,1);
lines.einsteinACoefficient = zeros(numberOfLines,1);
lines.airBroadenedWidth = zeros(numberOfLines,1);
lines.selfBroadenedWidth = zeros(numberOfLines,1);
lines.lowerStateEnergy = zeros(numberOfLines,1);
lines.temperatureDependence = zeros(numberOfLines,1);
lines.pressureShift = zeros(numberOfLines,1);
lines.upperVibrationalQuanta = char( zeros(numberOfLines, 15) );
lines.lowerVibrationalQuanta = char( zeros(numberOfLines, 15) );
lines.upperLocalQuanta = char( zeros(numberOfLines, 15) );
lines.lowerLocalQuanta = char( zeros(numberOfLines, 15) );
lines.errorCodes = char( zeros(numberOfLines, 6) );
lines.referenceCodes = char( zeros(numberOfLines, 12) );
lines.flagForLineMixing = char( zeros(numberOfLines, 1) );
lines.upperStatisticalWeight = zeros(numberOfLines, 1);
lines.lowerStatisticalWeight = zeros(numberOfLines, 1);

waitbar(0, h, 'Loading lines...');
tic % Start a timer.
for n=1:numberOfLines
    lines.moleculeNumber(n,1) = str2double( fread(fileID, 2, 'char=>char') );
    lines.isotopologueNumber(n,1) = str2double( fread(fileID, 1, 'char=>char') );
    lines.transitionWavenumber(n,1) = str2double( fread(fileID, 12, 'char=>char') );
    lines.lineIntensity(n,1) = str2double( fread(fileID, 10, 'char=>char') );
    lines.einsteinACoefficient(n,1) = str2double( fread(fileID, 10, 'char=>char') );
    lines.airBroadenedWidth(n,1) = str2double( fread(fileID, 5, 'char=>char') );
    lines.selfBroadenedWidth(n,1) = str2double( fread(fileID, 5, 'char=>char') );
    lines.lowerStateEnergy(n,1) = str2double( fread(fileID, 10, 'char=>char') );
    lines.temperatureDependence(n,1) = str2double( fread(fileID, 4, 'char=>char') );
    lines.pressureShift(n,1) = str2double( fread(fileID, 8, 'char=>char') );
    lines.upperVibrationalQuanta(n,:) = fread(fileID, 15, 'char=>char');
    lines.lowerVibrationalQuanta(n,:) = fread(fileID, 15, 'char=>char');
    lines.upperLocalQuanta(n,:) = fread(fileID, 15, 'char=>char');
    lines.lowerLocalQuanta(n,:) = fread(fileID, 15, 'char=>char');
    lines.errorCodes(n,:) = fread(fileID, 6, 'char=>char');
    lines.referenceCodes(n,:) = fread(fileID, 12, 'char=>char');
    lines.flagForLineMixing(n,:) = fread(fileID, 1, 'char=>char');
    lines.upperStatisticalWeight(n,1) = str2double( fread(fileID, 7, 'char=>char') );
    lines.lowerStatisticalWeight(n,1) = str2double( fread(fileID, 7, 'char=>char') );
    if fread(fileID,1) ~= 13
        error('The line endings of the HITRAN par file should be CR+LF, however the CR was missing.');
    end
    if fread(fileID,1) ~= 10
        error('The line endings of the HITRAN par file should be CR+LF, however the LF was missing.');
    end
    
    % Update the progress bar after every thousand lines have been read.
    if mod(n, 1000) == 0
        elapsedTimeInMinutes = toc / 60;
        linesRemaining = numberOfLines - n - 1;
        linesPerMinute = 1000 / elapsedTimeInMinutes;
        minutesRemaining = ceil( linesRemaining / linesPerMinute );
        waitbar(n / numberOfLines, h, sprintf('Loading lines...%d min remaining.', minutesRemaining) );
        tic % Reset the timer.
    end    
end

% Close the waitbar.
close(h);

% Close the file.
fclose(fileID);


