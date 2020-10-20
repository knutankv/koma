function [U_ord, sigma_ord, V_ord] = sort_fdd(U, sigma, V, phi_ref, min_mac)

[n_channels, n_modes] = size(phi_ref);
n_freqs = size(sigma, 3);

U_ord = nan(n_channels, n_modes, n_freqs);
V_ord = nan(n_channels, n_modes, n_freqs);
sigma_ord = repmat(diag(nan(n_modes, 1)), [1,1,n_freqs]);

for k = 1:size(sigma,3)
    mac = koma.modal.xmacmat(U(:, :, k), phi_ref);   
    [largest_macs, sv_ix_per_mode] = max(mac, [], 1);
    
    sv_ix = [];
    md_ix = [];
    sv_groups = unique(sv_ix_per_mode);

    for sv_group = sv_groups
        ix = find(sv_ix_per_mode==sv_group);
        [max_val, ixix] = max(largest_macs(ix));
        current_sv_ixs = sv_ix_per_mode(ix);
        if max_val > min_mac
            sv_ix = [sv_ix, current_sv_ixs(ixix)];
            md_ix = [md_ix, ix(ixix)];
        end
    end
    
    U_ord(:, md_ix, k) = U(:, sv_ix, k);
    V_ord(:, md_ix, k) = V(:, sv_ix, k);
    sigma_ord(md_ix, md_ix, k) = sigma(sv_ix, sv_ix, k);
end
    