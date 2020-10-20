function [phi] = maxreal(phi)
    for mode = 1:size(phi,2)
        angleVector=0:0.01:pi/2;
        rot_mode=phi(:,mode)*exp(angleVector*1i);
        [~,iReMax] = max(sum(real(rot_mode).^2,1));  %find index with the angle giving the largest real part
        phi(:,mode)=phi(:,mode)*exp(angleVector(iReMax)*1i);
        phi(:,mode) = phi(:,mode)*sign(sum(real(phi(:,mode))));     %make net positive
    end
end

