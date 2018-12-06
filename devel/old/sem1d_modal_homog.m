% [eigval,err,eigvec,x]=sem1d_homog(nel,ngll,bc,fig)
%
% Driver for sem1d_modal_analysis.m on a homogeneous 1D medium
% Solved on a regular mesh
%
% INPUT		nel	total number of elements
%		ngll	numer of GLL nodes per element
%		bc	'N' for Neumann boundary condition
%			'D' for Dirichlet
%		fig	[0] plots eigenfrequencies
%
% OUTPUT	eigval	eigenfrequencies, in ascending order
%		err	absolute error on eigenvalues ~ err/eigval
%		eigvec	eigenvectors (modes)
%		x	physical coordinates of nodes
%
function [eigval,err,eigvec,x]=sem1d_homog(nel,ngll,bc,fig)

if ~exist('fig','var'), fig=0; end

xmesh=(0:nel)' /nel;
rho=1;
mu=1;

if bc=='N'
  [eigvec,eigval,x,err]=sem1d_modal_analysis(xmesh,ngll,'NN',rho,mu);
  k=(0:length(eigval)-1)';
else
  [eigvec,eigval,x,err]=sem1d_modal_analysis(xmesh,ngll,'DD',rho,mu);
  k=(1:length(eigval))';
end
anaval=pi*k;

if fig, plot(k,eigval,'o-',k,anaval,'--'); end

eigval= [eigval,anaval];


%--------
% Personal notes
% 
% for k=3:10,
%   [val,vec,x]=sem1d_homog(1,k);
%   hold on
%   wmax(k)=max(val(:,1));
%   dx(k)=x(2);
% end
% loglog(dx,wmax,'o',[0.03:0.01:0.5],2.33./[0.03:0.01:0.5])
%
% ... so: wmax = 7/3 * c/dx
%         fmax = 7/3 /(2*pi) *c/dx
% also dx = 4/(ngll^2-1)
% See email to Matt Haney Mon Aug 23 18:41:14 EDT 2004
% For stability of a discrete time scheme, we must have:
%   dt*wmax < critical_CFL
% For non-dissipative Newmark-alpha critical_CFL=2, so:
%   dt*c/dx < 2* 3/7
%   CFL < 0.8571

