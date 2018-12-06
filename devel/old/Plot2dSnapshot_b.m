% Plot2dSnapshot_b(x,y,v,indx,vsat)
%
% INPUT		x,y	2D coordinates of GLL nodes (global vector)
%		v	field to be plotted (global vector)
%		indx	cell info, output from Init2dSnapshot
%		vsat	color scale [-vsat vsat]
% OUTPUT	h	patch handle 
%
function h=Plot2dSnapshot_b(x,y,v,indx,vsat)

[NGLL,NGLL,NELX,NELY] = size(x);
NEL=NELX*NELY;

clf
set(gcf,'DefaultPatchEdgeColor','none');
axis([x(1,1,1) x(NGLL,1,NELX) y(1,1,1) y(NGLL,NGLL,end)])
for e=1:NEL,
  xl = x(:,:,e);
  yl = y(:,:,e);
  vl = v(:,:,e);
  h=patch('Vertices',[xl(:) yl(:)],'Faces',indx, ...
          'FaceVertexCData',vl(:),'FaceColor','interp');
end
axis equal tight
title('SEM2D snapshot')
xlabel('X')
ylabel('Y')
if exist('vsat','var'), set(gca,'CLim',[-vsat vsat]); end
colorbar('vert')
