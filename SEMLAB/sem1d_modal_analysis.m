% Modal analysis of a 1D medium discretized by the Spectral Element Method
% Scalar wave equation
%
% [eigvec,eigval,x,err]=sem1d_modal_analysis(X,NGLL,BC,RHO,MU)
%
% INPUT		X(:)	physical coordinates of the elements vertices (ordered)
%		NGLL	number of Gauss-Lobatto-Legendre points per element 
%			= polynomial degree +1
%		BC	['NN'] boundary conditions at min(X) and max(X) respectively
%				N	homogeneous Neumann (stress free)
%				D	homogeneous Dirichlet (no displacement)
%		RHO	structure for the generation of density distribution
%			(see BuildKM_1d.m) 
%		MU	same for the shear modulus
%
% OUTPUT	eigvec	eigenvectors
%		eigvel	eigenvalues
%		x	coordinates of the GLL nodes
%		err	absolute error on eigenvalues ~ err/eigval
%
function [eigvec,eigval,x,err]=sem1d_modal_analysis(X,NGLL,BC,RHO,MU)

% Defaults: 
% homogeneous Neumann (stress free) on both ends
if ~exist('BC','var') || isempty(BC), BC='NN'; end 
% unit properties
if ~exist('RHO','var') || isempty(RHO), RHO=1; end 
if ~exist('MU','var') || isempty(MU), MU=1; end 

%---- Problem setup ----

[iglob,x,nglob] = mesh1d(X,NGLL);	
[M,Kloc] = BuildMK_1d(x,iglob,RHO,MU);

% Assemble the global stiffness matrix
% Note: for large scale problems it is smarter to store K as sparse
NEL = length(X)-1;	% total number of elements
K = zeros(nglob,nglob);
for e=1:NEL,
  ig = iglob(:,e);
  K(ig,ig) = K(ig,ig) + Kloc(:,:,e);
end 

%---- Boundary Conditions ----

% modify matrices if any Dirichlet condition: 
% eliminate the boundary degree of freedom
if ~strcmp('NN',BC)
  if BC(1)=='D', i1=2; else, i1=1; end
  if BC(2)=='D', i2=nglob-1; else, i1=nglob;  end
  M = M(i1:i2);
  K = K(i1:i2,i1:i2);
  nglob = i2-i1+1;
end


%---- Eigenvalue solver ----

% We're set up now to compute the eigenvalues and eigenvectors of M^(-1)*K 
% However, the eigenvalue solver works best with symmetric matrices
% and M^(-1)*K is NOT symmetric.
% The trick is to define a new symmetric matrix:   A = M^(-1/2)*K*M^(-1/2)
% The eigenvalues of the original problem are the same as those of A,
% the eigenvectors are V = M^(-1/2) * V_of_A

% building A = M^(-1/2)*K*M^(-1/2)
invsqrtM = 1./sqrt(M);
A = K .* repmat(invsqrtM',nglob,1);
A = repmat(invsqrtM,1,nglob) .*A;
A = 0.5*(A+A'); % just to be sure A is exactly symmetric

% Find all eigenvalues/eigenvectors.
% The standard Matlab routine 'eig' does not exploit sparseness :(
% The routine 'eigs' does, but only if a few eigenpairs are requested
[eigvec,eigval]=eig(A);

eigval=sqrt(abs(diag(eigval))); 	% eigenvalues of A are frequency^2

eigvec=repmat(invsqrtM,1,nglob) .* eigvec;	% rescale eigenvectors

% output in ascending eigenvalue order
[eigval,isort]=sort(eigval);
eigvec= eigvec(:,isort);

% if any Dirichlet condition: add the boundary node (=0) to eigenvectors
if ~strcmp('NN',BC)
  pad = zeros(1,nglob);
  if BC(1)=='D', eigvec = [ pad; eigvec ]; end 
  if BC(2)=='D', eigvec = [ eigvec; pad ]; end
end

% Very low frequency eigenfrequencies are not well recovered
% due to round-off error.
% An estimate of the absolute error on eigenvalues is err./eigval
if nargout>=4, err=0.5*eps*norm(A,2); end
