function [B] = findDampingBasis(xd)
    if length(xd)==2
         y1 = 1;
         y2 = -xd(1)/xd(2);
         y = [y1;y2];
         B = [xd./norm(xd), y./norm(y)];
    elseif length(xd)==3
        %%% here--- try to using the xd give a damping matrix
        %      using same algorithm as the HOIMDS --- Gram–Schmidt process
         Z = [0 0 1]';
         eig1=xd;
         middle1=Z'*eig1;
         middle2=eig1'*eig1;
         Z=Z-(middle1/middle2)*eig1;
         Z=Z / norm(Z);
         eig2=cross(eig1,Z);
         eig2=eig2 / norm(eig2);
         
         B = [eig1, eig2 , Z];
    else
        msg = 'findDampingBasis error.';
        error(msg)
    end
end
