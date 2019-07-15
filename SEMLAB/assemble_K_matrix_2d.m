function [ Kglob ] = assemble_K_matrix_2d(NelX, NelY, NGLL, dxe, dye, nglob, iglob,W)
%ASSEMBLE_K_MATRIX_2D Summary of this function goes here
%   Detailed explanation goes here

    Nel = NelX*NelY;
    
    %--------------------------------------
    % Element wise local stiffness matrix:
    %--------------------------------------
    
    [xgll,wgll,H] = GetGLL(NGLL);
    wgll2 = wgll * wgll' ;
    
    dx_dxi  = 0.5*dxe;
    dy_deta = 0.5*dye;
    jac = dx_dxi*dy_deta;
    
    w_temp = zeros(NGLL,NGLL); % material properties
    Ke_temp = zeros(NGLL,NGLL); % temporary matrix
    
    % Stiffness matrix for each element
    %Ke = zeros(NGLL*NGLL,NGLL*NGLL,Nel);
    Ke = zeros(NGLL*NGLL,NGLL*NGLL);
    
    % Matrix assembly
    I = zeros(length(iglob(:)),1);
    J = zeros(length(iglob(:)),1);
    V = zeros(length(iglob(:)),1);
    
    ct = 1;
    
    del = eye(NGLL); % identity matrix: kronecker delta
    
    for eo = 1:Nel
        Ke_temp(:,:) = 0.;
        w_temp = W(:,:,eo);
        
        for i = 1:NGLL
            for j = 1:NGLL
                for k = 1:NGLL
                    for l = 1:NGLL
                        term1 = 0.; term2=0.;
                        for p = 1:NGLL
                           term1 = term1 + del(i,k)*w_temp(k,p)*(jac/dy_deta^2)*H(j,p)*H(l,p);
                           term2 = term2 + del(j,l)*w_temp(p,j)*(jac/dx_dxi^2)*H(i,p)*H(k,p);
                        end
                        Ke_temp(i,j,k,l) = term1 + term2;
                    end
                end
            end
        end
        
        Ke(:,:) = reshape(Ke_temp, NGLL*NGLL,NGLL*NGLL);
            
        % Assembly
        % I,J,V are such that K(I,J) = V
        ig = iglob(:,:,eo);
        ig = ig(:);
        
        for j = 1:length(ig)
            for i = 1:length(ig)
                I(ct) = ig(i);
                J(ct) = ig(j);
                V(ct) = Ke(i,j);
                ct = ct + 1;
            end
        end
    end
    
    % Faster assembly as SparseCSC format
    Kglob = sparse(I,J,V,nglob,nglob);
    
    % Naive assembly: super slow for large systems
    % Kglob = spalloc(nglob,nglob,NGLL*nglob); % preallocate sparse matrix of size
                                               % nglob x nglob with space
                                               % for NGLL x nglob nonzeros
    % for eo = 1:Nel
    %    ig = iglob(:,:,eo);
    %    Kglob(ig(:),ig(:)) =  Kglob(ig(:), ig(:)) + Ke(:,:,eo);
    % end
    
end

