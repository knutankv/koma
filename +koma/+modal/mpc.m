function MPC = mpc(phi)

n_modes = size(phi, 2);
MPC = zeros(n_modes, 1);

for mode = 1:n_modes
    phin = phi(:, mode);
    Sxx = dot(real(phin), real(phin));
    Syy = dot(imag(phin), imag(phin));
    Sxy = dot(real(phin), imag(phin));
    
    eta = (Syy-Sxx)/(2*Sxy);
    
    lambda1 = (Sxx+Syy)/2 + Sxy*sqrt(eta^2+1);
    lambda2 = (Sxx+Syy)/2 - Sxy*sqrt(eta^2+1);
    
    MPC(mode) = ((lambda1-lambda2)/(lambda1+lambda2))^2;
end