function [values, has_qc] = read_arm_qc0(file, var_name)
% READ_ARM_QC0  Read an ARM variable and keep only its best-quality samples.
%
%   [values, has_qc] = read_arm_qc0(file, var_name)
%
%   Most ARM variables come with a quality-control companion "qc_<variable>".
%   Each bit of the QC value records one failed test (e.g. "value below the
%   valid minimum"); ARM rates each test "Bad" or "Indeterminate".
%   QC = 0 means the sample passed EVERY test. This function keeps only
%   those samples, the strictest choice (it also drops "Indeterminate"
%   samples); all others are set to NaN.
%
%   If the variable has no qc_ companion, nothing is removed and has_qc is
%   false (the variable is used as distributed).
%
%   Inputs   file      path to a netCDF file
%            var_name  variable name, e.g. 'downwelling_longwave'
%   Outputs  values    column vector of doubles (NaN = missing or failed QC)
%            has_qc    true if a qc_ field existed
%
%   Note: ncread already turns float fill values (_FillValue) into NaN.

values = double(ncread(file, var_name));
values = values(:);

info   = ncinfo(file);
names  = {info.Variables.Name};
has_qc = any(strcmp(names, ['qc_' var_name]));

if has_qc
    qc = double(ncread(file, ['qc_' var_name]));
    values(qc(:) ~= 0) = NaN;      % keep QC == 0 only
end

end
