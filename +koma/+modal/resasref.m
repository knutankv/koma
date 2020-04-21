function [lambda_save, phi_save, okrefidx, okididx] = resasref(lambda, phi, phi_ref, maccrit)
%% Restructure modes as reference (shuffle the mode numbering).
% 
% Arguments
% ---------------------------
% lambd : double
%     array of all eigenvalues
% phi : double
%     modal transformtaion matrix (all mode shapes), stacked column-wise
% phi_ref : int
%     modal transformation matrix from reference (the one to sort after)
% maccrit : int
%   MAC criterion (nan/empty if not satisfied)
%
% Returns
% ---------------------------
% lambda_save : double
%     output lambda (picked and sorted from lambda)
% phi_save : double
%     output phi (picked and sorted from phi)
% okrefidx : int
%     indices related to reference that are deemed ok
% okididx : int
%     indices related to input phi that are deemed ok


macs = modal.mac_numbers(phi_ref,phi);
[maxmacs,idxmodes]=max(macs);
okididx=find(maxmacs>maccrit);  % ok indices of identified
okrefidx = idxmodes(maxmacs>maccrit);   % ok indicies of reference
    
lambda_save = nan(1,size(phi_ref,2));
lambda_save(okrefidx) = lambda(okididx);
phi_save = nan(size(phi_ref));
phi_save(:,okrefidx) = phi(:,okididx);