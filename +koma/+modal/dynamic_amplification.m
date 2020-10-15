function D = dynamic_amplification(xi, varargin)

if nargin==1
    D = @(beta) (sqrt((1-beta.^2).^2+(2*xi*beta.^2).^2)).^(-1);
else
    omega0 = varargin{1};
    D = @(omega) (sqrt((1-(omega/omega0).^2).^2+(2*xi.*(omega./omega0).^2).^2)).^(-1);

end
