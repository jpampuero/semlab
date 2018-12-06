%GETGLL Gauss-Lobatto-Legendre points, weights and derivatives of Lagrange polynomials 
%
% [x,w,h] = GetGLL(n)
%
% INPUT		n	Number of GLL points per 1D spectral element 
%			= polynomial degree +1
%		
% OUTPUT	x(n)	Coordinates of GLL points in the reference segment [-1:1]
%		w(n)	Quadrature weights
%		h(n,n) 	Derivatives of Lagrange polynomials at the GLL nodes
%			h(i,j) = h'_i( x(j) )
%
% NOTE	Reads from pre-tabulated files for the usual range ( NGLL=3:20 )
% 	located in the directory "gll_xwh/" 
%
function [x,w,h]=GetGLL(ngll,kind)

if nargin>1, 
  prefix=kind(1:3);
else
  prefix = 'gll';
end
name = sprintf('gll_xwh/%s_%0.2u.tab',prefix,ngll);

% get full path
pathstr = fileparts(mfilename('fullpath'));
name = fullfile(pathstr,name);

if ~exist(name,'file')
  error(sprintf('Data file %s does not exist',name))
end

fid=fopen(name);
data=fscanf(fid,'%f',[ngll,ngll+2]);
fclose(fid);

x=data(:,1);
w=data(:,2);
h=data(:,3:end)';
