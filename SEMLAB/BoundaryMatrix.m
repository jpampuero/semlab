%BOUNDARYMATRIX builds the database for a boundary of the 2D SEM grid
%
% [B,iB,jB] = BoundaryMatrix(wgll,NELXY,iglob,jac1D,side)
%
% INPUT		wgll	GLL weights (see GetGLL)
%		NELXY(2) number of elements along each direction
%		iglob	local-to-global index table of the SEM grid 
%		jac1D	line jacobian
%		side	'L'	left
%			'R'	right
%			'T'	top
%			'B'	bottom
%
% OUTPUT	B	boundary matrix
%		iB	boundary-to-bulk global index table:
%			the k-th node of the boundary mesh
%			is the iB(k)-th node of the mesh
%		jB	boundary local-to-global index table:
%			the k-th GLL node of the e-th boundary element
%			is the jB(k,e)-th node of the boundary mesh
%
% NOTE	jB is only needed for interpolation on non-GLL nodes (see interp1_sem)
%
function [B,iB,jB] = BoundaryMatrix(wgll,NELXY,iglob,jac1D,side)

NGLL = length(wgll);
NELX = NELXY(1);
NELY = NELXY(2);

switch side

  case 'L'
    eB = (0:NELY-1)*NELX+1;
    igll = 1;
    jgll = 1:NGLL;

  case 'R'
    eB = (0:NELY-1)*NELX+NELX;
    igll = NGLL;
    jgll = 1:NGLL;

  case 'T'
    eB = (NELY-1)*NELX+ (1:NELX);
    igll = 1:NGLL;
    jgll = NGLL;

  case 'B'
    eB = 1:NELX;
    igll = 1:NGLL;
    jgll = 1;

end

NELB = length(eB);
ng = NELB*(NGLL-1)+1;
iB = zeros(ng, 1); 
B = zeros(ng, 1); 
jB = zeros(NGLL,NELB);
for e=1:NELB,
  ip = (NGLL-1)*(e-1)+[1:NGLL];
  iB(ip) = iglob(igll,jgll,eB(e));
  jB(:,e) = ip;
  B(ip) = B(ip) + jac1D*wgll;
end
