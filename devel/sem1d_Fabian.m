% SEM1D	applies the Spectral Element Method
% to solve the 1D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source.
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X
L=10;
NEL = 10;
X = [0:NEL]'*L/NEL;

%------------------------------------------
% STEP 2: INITIALIZATION

P = 6; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 1000; % number of timesteps

ABSO_BOTTOM =1;

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);

iglob = zeros(NGLL,NEL);	% local to global index mapping
rho = zeros(NGLL,NEL);		% density
mu = zeros(NGLL,NEL);		% shear modulus
nglob = NEL*(NGLL-1) + 1;	% number of global nodes
coor = zeros(nglob,1);		% coordinates of GLL nodes
M = zeros(nglob,1);		% global mass matrix, diagonal
K = zeros(NGLL,NGLL,NEL);	% local stiffness matrix
CFL = 0.85; 			% stability number
dt = Inf;  			% timestep (set later)

for e=1:NEL, % FOR EACH ELEMENT ...

 % The table I = iglob(i,e) maps the local numbering of the 
 % computational nodes (i,e) to their global numbering I.
 % 'iglob' is used to build global data 
 % from local data (contributions from each element)
  iglob(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';

 % Coordinates of the computational (GLL) nodes
  dxe = X(e+1)-X(e);
  coor(iglob(:,e)) = 0.5*(X(e)+X(e+1)) + 0.5*dxe*xgll ;

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
  rho(:,e) = 1;
  mu(:,e) = 1;

 % For this simple mesh the jacobian of the global-local
 % coordinate map is a constant
  dx_dxi = 0.5*dxe;

 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
  M(iglob(:,e)) = M(iglob(:,e)) + wgll .*rho(:,e) *dx_dxi;

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  W = mu(:,e).*wgll/dx_dxi;
%  K(:,:,e) = H * diag(W)* H';
  K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
  vs = sqrt(mu(:,e)./rho(:,e)); 
  vs = max( vs(1:NGLL-1), vs(2:NGLL) );
  dx = abs(diff( coor(iglob(:,e)) ));
  dt = min(dt, min(dx./vs));
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;

% absorbing boundaries
% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
temp1 = M/dt;
temp2 = temp1;
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  temp1(1) = temp1(1) + BcTopC/2;
  temp2(1) = temp2(1) - BcTopC/2;
end

% Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
va = zeros(nglob,1);
a = zeros(nglob,1);
f = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
F_IS_WAVE = 1;
if ~F_IS_WAVE,
  Fx = 5;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 0.25; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function, at t-dt
Ft = ricker( (0:NT-1)'*dt, Ff0,Ft0);

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe:L]';
OUTnx = length(OUTx);
OUTix = zeros(OUTnx,1);
% relocate to nearest GLL node
OUTdist = zeros(OUTnx,1);
for i=1:OUTnx,
  [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) );
end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,OUTnt);
OUTv = zeros(OUTnx,OUTnt);
OUTa = zeros(OUTnx,OUTnt);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Central difference scheme, as in NERA 
% http://gees.usc.edu/GEES/Software/NERA/2001/Default.htm

% va = mid-step velocity
% initial conditions: d=0 and va=0

for it=1:NT,

 % internal forces using previous displacement -K*d(t-1) 
  f(:) = 0; % store -K*d in a global array
  for e=1:NEL,
    ix = iglob(:,e);
    f(ix) = f(ix) - K(:,:,e)*d(ix) ;
  end 

 % add external forces
 % absorbing boundary and incoming wave
  if F_IS_WAVE 
    if ABSO_BOTTOM, 
      f(1) = f(1) +2*BcBottomC*Ft(it); 
    else
      f(1) = f(1) +BcBottomC*Ft(it); 
    end
  else
    f(Fix) = f(Fix) + Ft(it); 
  end

  v = va; % temporarily store va(t-1) in v

 % update va(t), from the viscoelastic dynamic equations
  va = ( va.*temp2 + f ) ./temp1 ; 

 % update d(t)
  d = d + dt*va;

 % note: at this point v contains va(t-1)
 % update a(t-1) and v(t-1) 
  a = (va - v)/dt;
  v = 0.5*(v + va);


%------------------------------------------
% STEP 4: OUTPUT
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;
    OUTd(:,OUTit) = d(OUTix);
    OUTv(:,OUTit) = v(OUTix);
    OUTa(:,OUTit) = a(OUTix);
% note: OUTa and OUTv are delayed by dt with respect to OUTd 
%       (see central difference scheme) 
  end

end % of time loop

t = (1:OUTit) *OUTdt*dt -dt; % -dt: see note above
PlotSeisTrace(OUTx,t,OUTv);
