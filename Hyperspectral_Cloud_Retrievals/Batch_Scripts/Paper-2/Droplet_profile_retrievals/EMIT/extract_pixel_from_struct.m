%% Extract a single pixel from a structure by dynamically finding pixel-dependent fields
%
% This function recursively walks through all fields of a structure using
% fieldnames(). For each numeric array, it checks whether any dimension
% matches total_pixels. If so, it extracts only the data for pixel pp along
% that dimension. If no dimension matches total_pixels, the field is kept
% as-is (it's global/constant data).
%
% This approach is future-proof: if new fields are added to the structure,
% they will be handled automatically without modifying this function.
%
% INPUTS:
%   s            - The structure to extract from
%   pp           - The pixel index to extract (1-based)
%   total_pixels - The total number of pixels in the trimmed data
%
% OUTPUTS:
%   s_out - A new structure with the same fields, where pixel-dependent
%           arrays have been reduced to a single pixel
%
% EXAMPLE:
%   emit_single = extract_pixel_from_struct(emit, 5, 30);
%
% By Andrew John Buggee
%%

function s_out = extract_pixel_from_struct(s, pp, total_pixels)


% Get all field names at this level
fields = fieldnames(s);

s_out = struct();

for ff = 1:length(fields)

    field_val = s.(fields{ff});

    if isstruct(field_val)

        % ----- Recurse into nested structs -----
        s_out.(fields{ff}) = extract_pixel_from_struct(field_val, pp, total_pixels);

    elseif isnumeric(field_val) || islogical(field_val)

        % ----- Numeric or logical array: check dimensions -----
        sz = size(field_val);

        % Find which dimension(s) match total_pixels
        matching_dims = find(sz == total_pixels);

        if ~isempty(matching_dims) && total_pixels > 1

            % Extract along the first matching dimension
            dim = matching_dims(1);

            % Build indexing expression: s_out.field = field_val(:, ..., pp, :, ...)
            idx = repmat({':'}, 1, ndims(field_val));
            idx{dim} = pp;
            s_out.(fields{ff}) = field_val(idx{:});

        else

            % No dimension matches total_pixels: keep as-is (global data)
            s_out.(fields{ff}) = field_val;

        end

    else

        % ----- Cell arrays, strings, datetime, etc.: keep as-is -----
        s_out.(fields{ff}) = field_val;

    end

end


end
