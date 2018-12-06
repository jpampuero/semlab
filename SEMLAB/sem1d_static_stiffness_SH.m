% Compute the static stiffness of an elastic medium with depth-dependent properties
% in anti-plane deformation
% as a function of horizontal wavenumber
% The stiffness is the amplitude ratio between a sinusoidal stress applied at the surface
% and the resulting surface displacement
% The bottom boundary has fixed displacement

function Kstatic = sem1d_static_stiffness_SH(k,X,NGLL,MU)

[iglob,x,nglob] = mesh1d(X,NGLL);	
[M,Kloc] = BuildMK_1d(x,iglob,MU,MU); 
K = AssembleMatrix1D(Kloc,iglob);

% enforce zero displacement at the bottom by removing the node from the unknowns
nglob = nglob-1;
M = diag(M(1:nglob));
K = K(1:nglob,1:nglob);

% a force applied at the surface
F = zeros(nglob,1);
F(1) = 1;

nk = length(k);
Kstatic = zeros(nk,1);
for ik = 1:nk,
  A = K + k(ik)^2 * M;
  d = A \ F;
  Kstatic(ik) = 1/d(1);
end
