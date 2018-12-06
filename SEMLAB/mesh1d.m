%MESH1D builds a Spectral Element grid for a 1D medium
% 	composed of elements with internal Gauss-Lobatto-Legendre (GLL) nodes.
%
% [iglob,coor] = mesh1d(X,NGLL)
%
% INPUT		X(NEL+1) sorted and non-redundant list of coordinates of 
%			the edge nodes of the elements in the "macro-mesh". 
%			Can be irregular.
%		NGLL    Number of GLL points per element 
%			= polynomial degree +1
% 
% OUTPUT	iglob(NGLL,NEL) an index table that
%			maps the local numbering of the GLL nodes (i,e)
%			into their global numbering I=iglob(i,e).
%			This table is required to assemble global data 
%			from local data (adding contributions from each element).
%		coor(:) Non-redundant list of coordinates of the global GLL nodes
%		nglob	Total number of global nodes
%
function [iglob,coor,nglob] = mesh1d(X,NGLL)

NEL = length(X)-1;		% number of elements
nglob = NEL*(NGLL-1) + 1; 	% number of global GLL nodes
iglob = zeros(NGLL,NEL);
coor = zeros(nglob,1);
xgll = GetGLL(NGLL); % reference GLL nodes, in [-1:1]

for e=1:NEL,
  iglob(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';
  dxe = X(e+1)-X(e);
  coor(iglob(:,e)) = 0.5*(X(e)+X(e+1)) + 0.5*dxe*xgll ;
end
