% Plot2dSnapshot(x,y,v,indx,vsat)
%
% INPUT		x,y	2D coordinates of GLL nodes (global vector)
%		v	field to be plotted (global vector)
%		indx	cell info, output from Init2dSnapshot
%		vsat	color scale [-vsat vsat]
% OUTPUT	h	patch handle 
%
function h=Plot2dSnapshot(x,y,v,indx,vsat)

clf
set(gca,'DefaultPatchEdgeColor','none');
h=patch('Vertices',[x(:) y(:)],'Faces',indx,'FaceVertexCData',v(:),'FaceColor','interp');
axis equal 
axis tight
title('SEM2D snapshot')
xlabel('X')
ylabel('Y')
if exist('vsat','var'), set(gca,'CLim',vsat); end
colorbar('vert')
