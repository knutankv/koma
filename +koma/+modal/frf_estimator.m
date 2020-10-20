function H = frf_estimator(S, in_channels, out_channels, estimator_type)

warning off
n_freqs = size(S,3);
H = zeros(length(out_channels), length(out_channels), n_freqs);

if strcmp(estimator_type, 'H1')
    for k = 1:size(S,3)
        H(:, :, k) = S(in_channels, in_channels, k)\S(in_channels, out_channels, k).';
    end
elseif strcmp(estimator_type, 'H2')
    for k = 1:size(S,3)
        H(:, :, k) = S(out_channels, in_channels, k)\S(out_channels, out_channels, k).';
    end
else
    error('Wrong type defined. Use "H1" or "H2"')
end

warning on