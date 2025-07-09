% Compute the information content gained by performing the retireval

function [H] = compute_shannonInformation_content(posterior_cov, prior_cov)


% compute the Shannon information content
H = -1/2 * log2(posterior_cov * prior_cov^(-1));


end