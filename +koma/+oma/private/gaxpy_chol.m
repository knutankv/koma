function G = gchol(A)
%% Cholesky factorization of a symmetric positive definite matrix. 
%
% Arguments
% -------------------
% A : double
%   n-by-n array which is symmetric and positive definite
% 
% Returns
% -------------------
% G : double
%   lower triangular Cholesky factor, A = GG'
%
% References
% -------------------
% Algorithm 4.2.1 from Golub et al. :cite:`Golub2012`

n = length(A);
v = zeros(n,1);

for j=1:n
    v(j:n) = A(j:n,j) - A(j:n,1:j-1)*A(j,1:j-1)';
    A(j:n,j) = v(j:n)/sqrt(v(j));
end

G = tril(A);