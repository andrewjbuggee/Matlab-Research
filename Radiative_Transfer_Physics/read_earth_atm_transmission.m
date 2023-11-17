%% Read Soruce file


%%

function [transmission] = read_earth_atm_transmission(filename)

fileID = fopen(filename);

format_spec = '%f %f';

Data = textscan(fileID, format_spec, 'Delimiter',' ',...
   'MultipleDelimsAsOne',1, 'CommentStyle','#');

transmission.wavelength = Data{1};
transmission.T = Data{2};



end

