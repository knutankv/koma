function M = smooth_3d(M, span)

for row = 1:size(M, 1)
    for col = 1:size(M, 2)
        M(row,col,:) = smooth(squeeze(M(row,col,:)), span);    
    end
end