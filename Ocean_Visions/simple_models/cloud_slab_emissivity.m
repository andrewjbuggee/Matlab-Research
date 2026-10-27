function [eps_cloud, eps_slab] = cloud_slab_emissivity(P)
% CLOUD_SLAB_EMISSIVITY  Emissivity of the cloud and of the atmospheric slab that
%                        holds it, from P.cloud.phase, P.cloud.WP_g_m2 and that
%                        slab's clear-air emissivity.
%
%   [eps_cloud, eps_slab] = cloud_slab_emissivity(P)
%
% The cloud particles and the layer's gases are independent absorbers, so their
% transmissions multiply:  eps_slab = 1 - (1 - eps_gas)(1 - eps_cloud).  Computed
% here, at run time, from the primary cloud variables -- so a caller that edits
% P.cloud.phase or P.cloud.WP_g_m2 between runs always gets the emissivity it
% asked for, whether or not it re-ran seb_column_finalize.
eps_cloud = cloud_lw_emissivity(P.cloud.WP_g_m2, P.cloud.phase, P);
if isfield(P.atm, 'cloud_layer') && ~isempty(P.atm.cloud_layer)
    eg = P.atm.layers(P.atm.cloud_layer).emis_gas;
    eps_slab = 1 - (1 - eg) * (1 - eps_cloud);
else
    eps_slab = NaN;
end
end
