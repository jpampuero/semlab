% [iglob,x,y]=MeshBox(LX,LY,NELX,NELY,XGLL)
% 
% PURPOSE	Generates a Spectral Element mesh for a rectangular box,
%		elements and their internal Gauss-Lobatto-Legendre (GLL) sub-grids.
%
% INPUT		LX	x-length of the box
%		LY	y-length of the box
%		NELX	number of elements along x
%		NELY	number of elements along y
%		NGLL	number of GLL nodes (polynomial degree +1)
%
% OUTPUT	iglob(NGLL,NGLL,NELX*NELY) maps the local numbering of the 
%			computational nodes to their global (non-redundant) numbering.
%			I = iglob(i,j,e) is the global node index of the (i,j)-th GLL
%			node internal to the e-th element
%			Elements are numbered row by row from bottom-left to top-right.
% 			The table iglob is tipically needed to build or assemble 
%			global data from local data (contributions from each element)
%		x(:)	global x-coordinates of the GLL nodes, start at 0
%		y(:)	global y-coordinates of the GLL nodes, start at 0
%
function [iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL)

dxe = LX/NELX;
dye = LY/NELY;

XGLL = GetGLL(NGLL);

iglob = zeros(NGLL,NGLL,NELX*NELY);	% local to global index mapping
nglob = (NELX*(NGLL-1)+1)*(NELY*(NGLL-1)+1);	% number of global nodes
x     = zeros(nglob,1);		% coordinates of GLL nodes
y     = zeros(nglob,1);	

e=0;
last_iglob = 0;
igL = reshape([1:NGLL*(NGLL-1)],NGLL-1,NGLL);
igB = reshape([1:NGLL*(NGLL-1)],NGLL,NGLL-1);
igLB = reshape([1:(NGLL-1)*(NGLL-1)],NGLL-1,NGLL-1);
xgll = repmat( 0.5*(1+XGLL) , 1,NGLL);
ygll = dye*xgll';
xgll = dxe*xgll;

for ey=1:NELY, 
for ex=1:NELX, 

  e = e+1;

 % Take care of redundant nodes at element edges :
  if e==1  % first element: bottom-left
    ig = reshape([1:NGLL*NGLL],NGLL,NGLL);
  else
    if ey==1 	%  elements on first (bottom) row
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(2:NGLL,:) = last_iglob + igL; 		% the rest
    elseif ex==1 % elements on first (left) column
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(:,2:NGLL) = last_iglob + igB; 		% the rest
    else 	% other elements
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(2:NGLL,2:NGLL) = last_iglob + igLB;
    end
  end
  iglob(:,:,e) = ig;
  last_iglob = ig(NGLL,NGLL);

 % Global coordinates of the computational (GLL) nodes
  x(ig) = dxe*(ex-1)+xgll;
  y(ig) = dye*(ey-1)+ygll;

end
end

