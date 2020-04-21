function [ macs ] = xmacmat(phi1, varargin)
%% Modal assurance criterion numbers, cross-matrix between two modal transformation matrices (modes stacked as columns).
% 
% Arguments
% ---------------------------
% phi1 : double
%     reference modes
% phi2 : double, optional
%     modes to compare with, if not given (i.e., equal default value None), phi1 vs phi1 is assumed
% conjugates : true, optional
%     check the complex conjugates of all modes as well (should normally be true)
% 
% Returns
% ---------------------------
% macs : boolean
%     matrix of MAC numbers


if nargin==1
    phi2 = phi1;
else
    phi2 = varargin{1};
end

norms = sum(phi1.'.*phi1',2)*sum(phi2.'.*phi2',2)';

if nargin==3 && varargin{2} == true
    macs = max(abs(phi1'*phi2).^2./norms,abs(phi1'*conj(phi2)).^2./norms);
else
    macs = abs(phi1'*phi2).^2./norms;
end

end

