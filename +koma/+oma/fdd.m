function [U, sigma, V, n_lines] = fdd(cpsd, n_lines)
[n_channels, ~, n_freqs] = size(cpsd);

if n_lines>n_channels
   disp(['Number of lines reduced to number of channels (', num2str(n_channels), ')']) 
   n_lines = n_channels;
end

U = zeros(n_channels, n_lines, n_freqs);
sigma = zeros(n_lines, n_lines, n_freqs);
V = zeros(n_channels, n_lines, n_freqs);

for k = 1:n_freqs
   [u, s, v] = svd(cpsd(:,:,k));
   U(:,:,k) = u(:, 1:n_lines);
   sigma(:,:,k) = s(1:n_lines, 1:n_lines);
   V(:,:,k) = v(:, 1:n_lines);   
end
