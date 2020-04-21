function [lambda_stab, phi_stab, order_stab, idx_stab] = find_stable_poles(lambda, phi, order, s, stabcrit, indicator)
%% Post-processing of Cov-SSI results, to establish modes (stable poles) from all poles.
% 
% Arguments
% ---------------------------
% lambda : double
%     cell array of arrays with complex-valued eigenvalues, one cell entry per order
% phi : double
%     cell array of matrices with complex-valued eigevectors (stacked column-wise), one cell entry per order
% orders : int
%     1d array with orders corresponding to the cell array entries in lambd and phi
% s : int
%     stability level, see :cite:`Kvale2017_OMA`
% stabcrit : {'freq': 0.05, 'damping':0.1, 'mac': 0.1}, optional
%     criteria to be fulfilled for pole to be deemed stable
% indicator : 'freq', optional
%     what modal indicator to use ('freq' or 'mac')
% 
% Returns
% ---------------------------
% lambd_stab : double
%     array with complex-valued eigenvalues deemed stable
% phi_stab : double
%    2d array with complex-valued eigenvectors (stacked as columns), each column corresponds to a mode
% order_stab : int
%     corresponding order for each stable mode
% idx_stab : int
%     indices of poles (within its order from the input) given in lambd_stab deemed stable (for indexing post)
% 
% References
% --------------------------------
% Kvaale et al. :cite:`Kvale2017_OMA`



wtol=stabcrit(1);
xitol=stabcrit(2);
mactol=stabcrit(3);

lambda_stab=[];
phi_stab=[];
order_stab=[];
idx_stab=[];

%% ESTABLISH FOR ORDER ABOVE STABLEVEL
for i = s+1:length(order)
    w = abs(lambda{i});
    xi = -real(lambda{i})./abs(lambda{i});
    phi0 = phi{i};    %modal transformation matrix for order with index i

    %STABLE POLES
    for m = 1:length(w)
        stab = 0;
        for level = 1:s
            philast = phi{i-level};

            switch indicator
                case 'mac'
                    macs = koma.modal.xmacmat(philast,phi0(:,m),true);
                    [~,mlast] = max(macs(:,1));   %find largest mac in first column (all in matrix compared with vector)
                case 'freq'
                    wlast=abs(lambda{i-level});
                    [~,mlast] = min(abs(w(m)-wlast));
            end

            lambdalast = lambda{i-level}(mlast);
            xilast = -real(lambdalast)./abs(lambdalast);
            wlast = abs(lambdalast);
            dxi = abs(xi(m)-xilast);
            dw = abs(w(m)-wlast);
            mac = koma.modal.xmacmat(philast(:,mlast),phi0(:,m),true);

            if (dw/wlast<=wtol)&&(dxi/xilast<=xitol) && (mac>=(1-mactol))
                stab=stab+1;
            else
                stab=0;
                break
            end
        end

        if  stab>=s
            lambda_stab = [lambda_stab; lambda{i}(m)];
            phi_stab =  [phi_stab, phi{i}(:,m)];
            order_stab = [order_stab; order(i)];
            idx_stab = [idx_stab; m];
        end
    end
end
