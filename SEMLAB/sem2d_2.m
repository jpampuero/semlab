% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid.
%
% Version 2: domain = rectangular
%            medium = general (heterogeneous)
%            boundaries = 1 fault + 3 paraxial
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
LX=30;
LY=10;
NELX = 60;
NELY = 20;
P = 8; % polynomial degree
%********

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
x = x-LX/2;
nglob = length(x);


%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
M     = zeros(nglob,1);		% global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

%**** Set here the parameters of the time solver : ****
NT = 2800; % number of timesteps
CFL   = 0.6; 			% stability number = CFL_1D / sqrt(2)
%********

dt    = Inf;  			% timestep (set later)

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
for ey=1:NELY, 
for ex=1:NELX, 

  e = (ey-1)*NELX+ex;
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  rho(:,:) = 1;
%  if ex>NELX/2-2 & ex<NELX/2+3
%    mu(:,:)  = 0.25;
%  else
    mu(:,:)  = 1;
%  end
%********

 % Diagonal mass matrix
  M(ig) = M(ig) + wgll2 .*rho *jac;

 % Local contributions to the stiffness matrix K
 %  WX(:,:,e) = wgll2 .* mu *jac/dx_dxi^2;
 %  WY(:,:,e) = wgll2 .* mu *jac/dy_deta^2;
  W(:,:,e) = wgll2 .* mu;

 % The timestep dt is set by the stability condition
 %   dt = CFL*min(dx/vs)
  vs = sqrt(mu./rho); 
  if dxe<dye
    vs = max( vs(1:NGLL-1,:), vs(2:NGLL,:) );
    dx = repmat( diff(xgll)*0.5*dxe ,1,NGLL); 
  else
    vs = max( vs(:,1:NGLL-1), vs(:,2:NGLL) );
    dx = repmat( diff(xgll)'*0.5*dye ,NGLL,1); 
  end
  dtloc = dx./vs;
  dt = min( [dt dtloc(1:end)] );

end
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;


%-- Absorbing boundaries (first order): 
[BcLC,iBcL] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'L'); % Left
[BcRC,iBcR] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'R'); % Right
[BcTC,iBcT] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi ,'T'); % Top
impedance = 1; % = sqrt(rho*mu)
BcLC = impedance*BcLC;
BcRC = impedance*BcRC;
BcTC = impedance*BcTC;

% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
% Fortunately C is diagonal.
M(iBcL) = M(iBcL) +half_dt*BcLC;
M(iBcR) = M(iBcR) +half_dt*BcRC;
M(iBcT) = M(iBcT) +half_dt*BcTC;


%-- DYNAMIC FAULT at bottom boundary
[FltB,iFlt] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi,'B');
FltN = length(iFlt);
FltZ = M(iFlt)./FltB /dt;
FltX = x(iFlt);
FltV = zeros(FltN,NT);
FltInitStress = repmat(0.55,FltN,1);
%FltInitStress( abs(FltX)<=2 ) = 0.601;
FltState = zeros(FltN,1);
FltFriction.MUs = repmat(0.6,FltN,1);
FltFriction.MUd = repmat(0.5,FltN,1);
FltFriction.MUs(abs(FltX)>10) = 100; % barrier
FltFriction.MUd(abs(FltX)>10) = 100; % barrier
FltFriction.Dc  = 0.1;
FltFriction.W = (FltFriction.MUs-FltFriction.MUd)./FltFriction.Dc;
FltStrength = friction(FltState,FltFriction) - FltInitStress; % strength excess


%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = [-5:0.25:5]';		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(4,OUTnseis,1);	% y coord of receivers
%********

[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 40;
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);



%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%
for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 

 % internal forces at mid-step -K*d(t+1/2) 
 % stored in global array 'a'
  a(:) = 0; 
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);
   %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
   %local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
    wloc = W(:,:,e);
    d_xi = wloc.*d_xi;
    d_xi = H * d_xi;
    d_eta = wloc.*d_eta;
    d_eta = d_eta *Ht;
    local = coefint1* d_xi  + coefint2* d_eta ;
   %assemble into global vector
    a(ig) = a(ig) -local;
  end 

 % absorbing boundaries:
  a(iBcL) = a(iBcL) - BcLC .* v(iBcL);
  a(iBcR) = a(iBcR) - BcRC .* v(iBcR);
  a(iBcT) = a(iBcT) - BcTC .* v(iBcT) ;

 % fault boundary condition: slip weakening
  FltVFree = v(iFlt) + dt*a(iFlt)./M(iFlt);
  TauStick = FltZ .*FltVFree;
  Tau = min(TauStick,FltStrength); 
  a(iFlt) = a(iFlt) - FltB .*Tau;

 % solve for a_new:
  a = a ./M ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;

  FltState = max(2*d(iFlt),FltState);
  FltStrength = friction(FltState,FltFriction)-FltInitStress;
% rupture nucleation through time weakening, like Andrews 76
  VNUC = 0.4;
%  if time(it)<10/VNUC
    %FltStrength = min( FltStrength,...
    ix = find(abs(FltX)<=10);
    FltStrength(ix) = min(FltStrength(ix), ...
              max(0.5,0.55+(abs(FltX(ix))-VNUC*time(it))*0.05/(1.0*dxe))...
                     - FltInitStress(ix) ) ;
%  end
  FltV(:,it) = 2*v(iFlt);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[0 0.1]);
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off

    drawnow

  end

end % ... of time loop

figure(3)
sdt=10;
sdx=4;
surf(time(1:sdt:end),FltX(1:sdx:end),FltV(1:sdx:end,1:sdt:end))
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')
