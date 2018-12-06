% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source,
% in a structured undeformed grid.
%
% Version 0: domain = square box (aspect ratio=1)
%            medium = homogeneous 
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

% This is usually done in two steps:
% the interval [0,LX]*[0,LY] is divided into NELX*NELY quadrangular elements
% (that's basically a quadrangular finite element mesh),
% then each element is sub-grided with (P+1)^2 Gauss-Lobatto-Legendre points 
% (GLL nodes) where P is the polynomial order of the spectral elements.
% Actually, in this example the mesh is so simple 
% that we do it in a single step without storing intermediate data.

%**** Set here the parameters of the square box domain and mesh : ****
LX=30;	% side-length
NELX = 30; % number of elements along each side
P = 6; % polynomial degree (inside each element, along each direction)
%********

LY=LX;
NELY = NELX;
dxe = LX/NELX;
dye = dxe;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element

% Generate the Spectral Element mesh
% The domain is partitioned into elements,
% each element contains a cartesian GLL subgrid
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x);

% The global numbering of the elements follows this convention:
%
%      ... ... ... ... NELX*NELY
% ^    ... ... ... ... ...
% | NELX+1 ... ... ... 2*NELX
% |     1   2  ... ... NELX  
% --->

% The local numbering of GLL nodes follows this convention:
% 
%  	(1,NGLL)...	(NGLL,NGLL)
% ^   	...	...	... 
% |	(1,1)	...	(NGLL,1)
% --->


%------------------------------------------
% STEP 2: INITIALIZATION

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL (3 to 10).
% xgll = location of the GLL nodes inside the reference segment [-1:1]
[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

%**** Set here the parameters of the time solver : ****
TMAX = 35; % total time
CFL   = 0.6; %0.6 % stability number = CFL_1D / sqrt(2)
%********

%**** Set here the physical properties of the homogeneous medium : ****
rho = 1;
mu  = 1;
vs = sqrt(mu/rho); 
%********

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
dxmin = x(2)-x(1);
dt = CFL*dxmin/vs;
NT = ceil(TMAX/dt); % number of timesteps

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant 
jac = (0.5*dxe)^2;
% dx_dxi = dy_deta
% coefint = jac/dx_dxi^2 = jac/dy_deta^2 = 1

% The diagonal mass matrix is stored in a global array. 
% It is assembled here, from its local contributions
% Nodes at the boundary between two or more elements get 
% contributions from all of them.
M = zeros(nglob,1);
for e=1:NEL, 
  ig = iglob(:,:,e);
  M(ig) = M(ig) + wgll2 *rho *jac;
end

% The stiffness matrix K is not pre-assembled
% We only store coefficients needed for its local contributions
W = wgll2 .* mu;

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;

%-- SOURCE TERM: point force, time function = Ricker wavelet
%**** Set here the source location : ****
Fx = 0; Fy = 0;
%********
[Fx,Fy,Fig] = FindNearestNode(Fx,Fy,x,y);
Ff0 = 0.5; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = src_timef( time-0.5*dt,'ricker', Ff0,Ft0);

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
theta=linspace(0,pi/2,21);
dist=20;
OUTxseis = dist*cos(theta); 
OUTyseis = dist*sin(theta); 
%OUTxseis = [0:dxe:LX]';		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
%OUTyseis = repmat(0,OUTnseis,1);	% y coord of receivers
%********
% receivers are relocated to the nearest node
% OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
theta = acos(OUTxseis/dist);
OUTv = zeros(OUTnseis,NT);
OUTd = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 50;
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);


%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%
half_dt = 0.5*dt;

for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 

 % internal forces at mid-step -K*d(t+1/2): 
  a(:) = 0; % store -K*d in a global array
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);
   %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
    d_xi  = W.*d_xi;
    d_eta = W.*d_eta;
    local = H*d_xi + d_eta*Ht ;
   %assemble into global vector
    a(ig) = a(ig) -local;
  end 

 % add external forces
  a(Fig) = a(Fig) + Ft(it);

 % acceleration: a = (-K*d +F)/M
  a = a ./M ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;


%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  OUTd(:,it) = d(OUTiglob);
  
  if mod(it,OUTdt) == 0 | it==NT
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[-0.5 0.5]);
    hold on
    plot(OUTxseis,OUTyseis,'^',Fx,Fy,'*')
    hold off

    drawnow

  end

end % ... of time loop
