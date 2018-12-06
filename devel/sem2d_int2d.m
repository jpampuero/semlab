% Compute the 2D spatial integral of some quantity over the whole domain
% Only for box domains meshed with a regular grid
%
% SYNTAX	If field is in local storage format:
%		  intg = sem2d_int2d(field,h)
%		If field is in global sorage format:
%  		  intg = sem2d_int2d(field,h,ibool)
%
% INPUT		field	in "local" format, size = [ngll,ngll,nelem]
%			in "global" storage format, size = [npoin]
%		h	size of spectral elements (assumed square)
%		ibool	local to global numbering, size=[ngll,ngll,nelem]
%			see SEM2D_READ_SPECGRID
%
% OUTPUT	intg	integral of field over the whole domain
%
function intg = sem2d_int2d(field,h,ibool)

if ndims(field)==3
  [ngll,ngll,nelem]=size(field);
else
  if ~exist('ibool','var'), error('sem2d_integral: if field is in global storage, ibool must be provided'); end
  [ngll,ngll,nelem]=size(ibool);
%else
%  error('sem2d_integral: invalid field size')
end
[xgll,wgll,H] = GetGLL(ngll);
wgll2 = wgll * wgll' ;

intg=0.0;
if ndims(field)==3
  for e=1:nelem,
    tmp = wgll2.*field(:,:,e);
    intg = intg + sum(tmp(:));
  end
else
  for e=1:nelem,
    k = ibool(:,:,e);
    tmp = wgll2.*field(k);
    intg = intg + sum(tmp(:));
  end
end

% For a box domain with regular mesh 
% the jacobian of global-local coordinate transformation is constant
jac = (0.5*h)^2;
intg = intg*jac;
