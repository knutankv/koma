function [S] = cpsd_from_fdd_result(U, sigma, V)

[n_channels, ~, n_freqs] = size(U);
S = zeros(n_channels, n_channels, n_freqs);

for k = 1:n_freqs
    S(:, :, k) = U(:, :, k) * sigma(:, :, k) * V(:, :, k)';
end