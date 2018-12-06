%BUILDMK_1D creates M and K matrices for 1D SEM
% 
% [M,K] = BuildMK_1d(coor,iglob [,rho,mu] );
%
% INPUT		coor	Non-redundant list of coordinates of the global GLL nodes
%		iglob	Local-to-global index map (see mesh1d)
%		rho,mu	Density and shear modulus 
%			Three input modes:
%			+ a constant
%			+ a matrix, in local storage (#igll,#elem)
%			+ a structure with fields:
%				data	all quantities needed to generate the values
%				fun	a function that generates the values
%					with the following syntax:
%					value = fun(data,e,x)
%					where e is an element index
%					and x(:) a list of node coordinates
%
% OUTPUT	M(NGLOB) Global mass matrix, diagonal 
%		K(NGLL,NGLL,NEL) Local (element-wise) stifness matrices
%
function [M,K] = BuildMK_1d(coor,iglob,rho,mu);

[NGLL,NEL] = size(iglob);
nglob = length(coor);

if nargin<3, rho=1; end
if nargin<4, mu =1; end

M = zeros(nglob,1);		% global mass matrix, diagonal
K = zeros(NGLL,NGLL,NEL);	% local stiffness matrix

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
[xgll,wgll,H] = GetGLL(NGLL);

for e=1:NEL, % for each element ...

  x = coor(iglob(:,e));

 % jacobian of the global-local coordinate map
  dx_dxi = 0.5*( x(NGLL)-x(1) );

 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
  rhol = GetLocalValue(rho,e,x);
  M(iglob(:,e)) = M(iglob(:,e)) + wgll .*rhol *dx_dxi;

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  mul = GetLocalValue(mu,e,x);
  W = mul.*wgll/dx_dxi;
  %K(:,:,e) = H * diag(W) * H';
  K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

end
