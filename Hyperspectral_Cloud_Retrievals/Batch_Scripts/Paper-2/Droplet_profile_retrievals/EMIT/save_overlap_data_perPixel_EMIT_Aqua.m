%% Save overlapping EMIT-Aqua pixel data as individual .mat files for parallel retrieval
%
% This function finds overlapping pixels between EMIT, MODIS, AIRS, and
% AMSR using findOverlap_pixels_EMIT_Aqua_coincident_data, trims
% unnecessary data, and saves one .mat file per EMIT pixel. Each .mat file
% contains all the information needed to run the droplet profile retrieval
% for a single EMIT pixel independently.
%
% This is designed to be run once per coincident data subdirectory (e.g.,
% '2023_9_16_T191118_1/'). The resulting .mat files can then be distributed
% across a SLURM job array where each task processes one pixel.
%
% INPUTS:
%   folder_paths - Structure containing paths (from define_EMIT_dataPath_and_saveFolders)
%       Must include:
%       .coincident_dataPath   - Parent path to coincident data
%       .coincident_dataFolder - Subdirectory name (with trailing slash)
%
%   print_status_updates - (logical) Print progress to console
%
% OUTPUTS:
%   n_pixels_saved - Number of .mat files saved
%   output_dir     - Path to the directory containing saved .mat files
%
% EXAMPLE USAGE:
%   folder_paths = define_EMIT_dataPath_and_saveFolders(0);
%   folder_paths.coincident_dataPath = '/path/to/coincident_EMIT_Aqua_data/southEast_pacific/';
%   folder_paths.coincident_dataFolder = '2023_9_16_T191118_1/';
%   [n_saved, out_dir] = save_overlap_data_perPixel_EMIT_Aqua(folder_paths, true);
%
% By Andrew John Buggee
%%

function [n_pixels_saved, output_dir] = save_overlap_data_perPixel_EMIT_Aqua(criteria, folder_paths, print_status_updates,...
    emit_pixels_per_modis, sub_directory)







%% Find overlapping pixels between EMIT and Aqua instruments

if print_status_updates == true
    disp([newline, 'Finding overlapping pixels for subdirectory: ',...
        sub_directory, newline])
end

folder_paths.coincident_dataFolder = sub_directory;

[overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
    criteria, emit_pixels_per_modis);

% If no pixels found, increase the horizontal inhomogeneity index
while isempty(overlap_pixels.modis.linear_idx) == true

    disp([newline, 'No overlapping pixels that meet defined criteria. Increasing H index....', newline])
    criteria.H = criteria.H + 0.1;

    [overlap_pixels, emit, modis, airs, amsr, folder_paths] = findOverlap_pixels_EMIT_Aqua_coincident_data(folder_paths,...
        criteria, emit_pixels_per_modis);

    if isempty(overlap_pixels.modis.linear_idx) == false

        disp([newline, 'Horizontal Inhomogeneity Index - H = ', num2str(criteria.H), newline])

    end

end


%% Remove data that is not needed

emit = remove_unwanted_emit_data(emit, overlap_pixels.emit);

modis = remove_unwanted_modis_data(modis, overlap_pixels.modis);

airs = remove_unwanted_airs_data(airs, overlap_pixels.airs);

amsr = remove_unwanted_amsr_data(amsr, overlap_pixels.amsr);


%% Create the output directory for per-pixel .mat files

% Remove trailing slash from coincident_dataFolder for use in filenames
subdir_name = sub_directory;
if subdir_name(end) == '/'
    subdir_name = subdir_name(1:end-1);
end

output_dir = [folder_paths.coincident_dataPath, 'overlap_data_per_pixel/'];

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end


%% Save one .mat file per EMIT pixel
%
% For each EMIT pixel pp, we must extract:
%   - EMIT data: pixel pp from the trimmed EMIT arrays (total_pixels entries)
%   - MODIS data: the ONE unique MODIS pixel that overlaps with EMIT pixel pp
%   - AIRS data:  the ONE unique AIRS pixel closest to EMIT pixel pp
%   - AMSR data:  the ONE unique AMSR pixel closest to EMIT pixel pp
%
% IMPORTANT: After remove_unwanted_*_data, the trimmed arrays have
% different lengths for each instrument:
%   - EMIT:  total_pixels entries (one per EMIT pixel)
%   - MODIS: n_unique_modis entries (one per unique MODIS pixel)
%   - AIRS:  n_unique_airs entries  (one per unique AIRS pixel)
%   - AMSR:  n_unique_amsr entries  (one per unique AMSR pixel)
%
% The overlap_pixels structure maps EMIT pixel pp to its corresponding
% MODIS/AIRS/AMSR pixel. We use the same unique/find logic that the
% retrieval functions use to determine which index in the trimmed array
% corresponds to each EMIT pixel.

total_pixels = length(overlap_pixels.emit.linear_idx);
n_pixels_saved = 0;

% --- Build the mapping from EMIT pixel index to unique MODIS/AIRS/AMSR index ---
% This replicates the logic in retrieve_dropProf_acpw_EMIT_Aqua_singlePix_ver2.m
% and write_AIRS_radiosonde_DAT_with_multiPixels.m

% MODIS: find the index of each EMIT pixel's MODIS pixel in the unique list
unique_modis_pix = unique(overlap_pixels.modis.linear_idx);
n_unique_modis = length(unique_modis_pix);
modis_idx_for_emit = zeros(total_pixels, 1);
for xx = 1:total_pixels
    modis_idx_for_emit(xx) = find(unique_modis_pix == overlap_pixels.modis.linear_idx(xx));
end

% AIRS: find the index of each EMIT pixel's AIRS pixel in the unique list
unique_airs_pix = unique(overlap_pixels.airs.linear_idx);
n_unique_airs = length(unique_airs_pix);
airs_idx_for_emit = zeros(total_pixels, 1);
for xx = 1:total_pixels
    airs_idx_for_emit(xx) = find(unique_airs_pix == overlap_pixels.airs.linear_idx(xx));
end

% AMSR: find the index of each EMIT pixel's AMSR pixel in the unique list
unique_amsr_pix = unique(overlap_pixels.amsr.linear_idx);
n_unique_amsr = length(unique_amsr_pix);
amsr_idx_for_emit = zeros(total_pixels, 1);
for xx = 1:total_pixels
    amsr_idx_for_emit(xx) = find(unique_amsr_pix == overlap_pixels.amsr.linear_idx(xx));
end


for pp = 1:total_pixels

    % Check if the EMIT pixel is masked out (all NaN)
    if all(isnan(emit.radiance.measurements(:, pp)))

        if print_status_updates == true
            disp([newline, 'EMIT pixel ', num2str(pp), ' is masked out. Skipping...', newline])
        end

        continue

    end

    % ---- Extract only pixel pp from EMIT ----
    % EMIT arrays have total_pixels entries after trimming.
    % extract_pixel_from_struct slices dimension matching total_pixels.
    emit_pp = extract_pixel_from_struct(emit, pp, total_pixels);                      %#ok<NASGU>

    % ---- Extract the one overlapping MODIS pixel ----
    % MODIS arrays have n_unique_modis entries after trimming.
    % modis_idx_for_emit(pp) gives the index into those trimmed arrays.
    modis_pp = extract_pixel_from_struct(modis, modis_idx_for_emit(pp), n_unique_modis); %#ok<NASGU>

    % ---- Extract the one overlapping AIRS pixel ----
    % AIRS arrays have n_unique_airs entries after trimming.
    airs_pp = extract_pixel_from_struct(airs, airs_idx_for_emit(pp), n_unique_airs);     %#ok<NASGU>

    % ---- Extract the one overlapping AMSR pixel ----
    % AMSR arrays have n_unique_amsr entries after trimming.
    amsr_pp = extract_pixel_from_struct(amsr, amsr_idx_for_emit(pp), n_unique_amsr);     %#ok<NASGU>

    % ---- Extract this pixel from overlap_pixels ----
    % overlap_pixels arrays all have total_pixels entries.
    overlap_pixels_pp = extract_pixel_from_struct(overlap_pixels, pp, total_pixels);     %#ok<NASGU>

    % Store the pixel number (always 1 in the saved file, since the
    % per-pixel structures contain only this one pixel's data)
    pixel_num = 1; %#ok<NASGU>

    % Define the output filename
    output_filename = [output_dir, 'overlap_EMIT_pixel_', num2str(pp, '%03d'), '_', subdir_name, '.mat'];

    % Save all data needed for the single-pixel retrieval
    save(output_filename, 'emit_pp', 'modis_pp', 'airs_pp', 'amsr_pp',...
        'overlap_pixels_pp', 'folder_paths', 'pixel_num', 'criteria', '-v7.3');

    n_pixels_saved = n_pixels_saved + 1;

    if print_status_updates == true
        disp(['Saved pixel ', num2str(pp), ' of ', num2str(total_pixels), ': ', output_filename])
    end

end


%% Print summary

if print_status_updates == true

    disp([newline, '--- Save Summary ---'])
    disp(['Subdirectory: ', sub_directory])
    disp(['Total EMIT pixels found: ', num2str(total_pixels)])
    disp(['Pixels saved (non-masked): ', num2str(n_pixels_saved)])
    disp(['Pixels skipped (masked): ', num2str(total_pixels - n_pixels_saved)])
    disp(['Output directory: ', output_dir])
    disp([newline])

end


end
