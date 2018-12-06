% PLOT2DSNAPSHOT plots a 2D SEM field (in global storage)
%
% indx = Plot2dSnapshot(iglob)		initialize
% h = Plot2dSnapshot(x,y,v,indx,vsat)	plot
%
% INPUT		iglob	local to global numbering table
%		x,y	2D coordinates of GLL nodes (global vector)
%		v	field to be plotted (global vector)
%		indx	cell info created by an "initialize" call
%		vsat	color scale [-vsat vsat]
%
% OUTPUT	indx	vertices of each GLL cell
%		h	patch handle 
%
function out = Plot2dSnapshot(x_iglob,y,v,indx,vsat)

if ~exist('vsat','var'), vsat=[]; end 
if nargin>1
  out = do_plot2d(x_iglob,y,v,indx,vsat);
else
  out = init_plot2d(x_iglob);
end

%------------ initialize -----------------
function indx = init_plot2d(iglob)

[NGLL,NGLL,NEL] = size(iglob);
indx = zeros(NEL*(NGLL-1)^2, 4);
ip=0;
for e=1:NEL,
  for i=1:NGLL-1,
  for j=1:NGLL-1,
    ip = ip+1;
    indx(ip,:) = [iglob(i,j,e) iglob(i+1,j,e) iglob(i+1,j+1,e) iglob(i,j+1,e)];
  end
  end
end

%--------------- plot -----------------------
function h = do_plot2d(x,y,v,indx,vsat)

clf
set(gca,'DefaultPatchEdgeColor','none');
h=patch('Vertices',[x(:) y(:)],'Faces',indx,'FaceVertexCData',v(:),'FaceColor','interp');
axis equal 
axis tight
title('SEM2D snapshot')
xlabel('X')
ylabel('Y')
if exist('vsat','var') && ~isempty(vsat), set(gca,'CLim',vsat); end
colorbar('vert')
