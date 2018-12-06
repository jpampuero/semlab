% Compute the static stiffness of an elastic medium with depth-dependent properties
% in anti-plane deformation
% as a function of horizontal wavenumber
% The stiffness is the amplitude ratio between a sinusoidal stress applied at the surface
% and the resulting surface displacement
% The bottom boundary has fixed displacement

function [Kstatic,ds] = sem1d_static_stiffness_PSV(ktab,X,NGLL,LAMBDA,MU)

[iglob,x,nglob] = mesh1d(X,NGLL);

% generate mu and lambda as (NGLL,NEL) tables
% ...

K10a = BuildKdd_1d('10',x,iglob,lambda); 
K00a = BuildKdd_1d('00',x,iglob,lambda+2*mu); 
K11a = BuildKdd_1d('11',x,iglob,mu); 
K01a = BuildKdd_1d('01',x,iglob,-mu); 

K10b = BuildKdd_1d('10',x,iglob,mu); 
K00b = BuildKdd_1d('00',x,iglob,-mu); 
K11b = BuildKdd_1d('11',x,iglob,-(lambda+2*mu)); 
K01b = BuildKdd_1d('01',x,iglob,-lambda); 

% enforce rigid bottom by removing the bottom node from the unknowns
ndof = nglob-1;

% a normal stress applied at the surface
F = zeros(ndof,1);
F(1) = 1;

nk = length(ktab);
Kstatic = zeros(nk,1);

if nargout>1, ds=zeros(ndof,nk); end

for ik = 1:nk,

  k = ktab(ik);
 
  % build the matrix
  A11loc = k^2*K00a + K11a;
  A21loc = k*(K10a + K01a);
  A12loc = k*(K10b + K01b);
  A22loc = k^2*K00b + K11b;
  A11 = AssembleMatrix1D(A11loc,iglob);
  A11 = A11(1:ndof,1:ndof);
  A21 = AssembleMatrix1D(A21loc,iglob);
  A21 = A21(1:ndof,1:ndof);
  A12 = AssembleMatrix1D(A12loc,iglob);
  A12 = A12(1:ndof,1:ndof);
  A22 = AssembleMatrix1D(A22loc,iglob);
  A22 = A22(1:ndof,1:ndof);
  A = [A11 A12; A21 A22];

  % solve
  d = A \ F;

  % store result
  Kstatic(ik) = 1/d(1);
  if nargout>1, ds(:,ik)=d; end

end
