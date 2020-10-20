function [U_ord, sigma_ord, V_ord] = sort_fdd(U, sigma, V, phi_ref)
U_ord = U*0;
V_ord = V*0;
sigma_ord = sigma*0;

for k = 1:size(sigma,3)
   mac = koma.modal.xmacmat(U(:,:,k), phi_ref);   
   [~, ix] = max(mac);
   
   U_ord(:, :, k) = U(:, ix, k);
   V_ord(:, :, k) = V(:, ix, k);
   sigma_ord(:, :, k) = sigma(ix, ix, k);
end
    