function phi = unique_modes(U, mac_max)
   
phi = [];

for k = 1:size(U,3)
    if size(phi, 2)>0
        macs = koma.modal.xmacmat(phi, U(:, 1, k));
    else
        macs = 0;
    end
    
    if max(macs(:))<mac_max
        phi = [phi, U(:, 1, k)];
    end
end
