function [ W ] = WMatrix(NelX, NelY, NGLL, wgll2, dxe, dye, rho, vs)
%WMATRIX Summary of this function goes here
%   Local contributions of stiffness matrix (mu*wgll2): the material
%   properties and the geometry will go here

    mu    = zeros(NGLL,NGLL);
    
    for ey = 1:NelY
        for ex = 1:NelX
            eo = (ey-1)*NelX + ex;
            
            % add here the properties of heterogeneous medium
            mu(:,:) = rho*vs^2;
            
            W(:,:,eo) = wgll2.*mu;
        end
    end
 
end

