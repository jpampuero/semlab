%FindNearestNode finds the nearest mesh node to a requested location (2D)
%
% [xout,yout,iglob,dist] = FindNearestNode(xin,yin,X,Y)
%
% INPUT		xin,yin	requested receiver locations, can be vectors
%		X,Y	global mesh node coordinates
%
% OUTPUT	xout,yout	relocated receiver locations	
%		iglob	global node indices of relocated receivers
%		dist	distance between requested and relocated receivers
%
function [xout,yout,iglob,dist] = FindNearestNode(xin,yin,X,Y)

nseis = length(xin);
dist = zeros(nseis,1);
iglob = zeros(nseis,1);
for k=1:nseis,
  [dist(k),iglob(k)] = min( (X-xin(k)).^2+(Y-yin(k)).^2 );
end
dist = sqrt(dist);
xout = X(iglob); 
yout = Y(iglob); 
