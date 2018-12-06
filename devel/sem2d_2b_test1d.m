% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid.
%
% Version 2b: domain = rectangular
%             medium = general (heterogeneous)
%             boundaries = 1 fault + 3 paraxial
%	      time scheme = explicit Newmark-alpha (with dissipation)
%
% THIS VERSION: compare to TestFlt1D of SEM2DPACK
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
LX=1;
LY=10;
NELX = 2;
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
NT = 800; % number of timesteps
CFL   = 0.5; % stability number = CFL_1D / sqrt(2)
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

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  rho(:,:) = 1;
%  if ex>NELX/2-2 & ex<NELX/2+3
%    mu(:,:)  = 0.25;
%  else
    mu(:,:)  = 1;
%  end

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
% Explicit Newmark-alpha scheme with dissipation
% alpha=1/2, gamma=1
r=0.5; % dissipation factor, in [0.5:1]
beta = 0.5 -2*r^2*(r-1)/(1+r)^3;
timco1 = dt^2*(1-2*beta)/2;
timco4 = beta*dt^2;

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;


%-- Absorbing boundaries (first order): 
impedance = 1; % = sqrt(rho*mu)
%% Left
%ng = NELY*(NGLL-1)+1;
%BcLeftIglob = zeros(ng,1);
%BcLeftC = zeros(ng,1);
%for ey=1:NELY,
%  ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
%  e=(ey-1)*NELX+1;
%  BcLeftIglob(ip) = iglob(1,1:NGLL,e);
%  BcLeftC(ip) = BcLeftC(ip) + dy_deta*wgll*impedance ;
%end
%% The mass matrix needs to be modified at the boundary
%% for the implicit treatment of the term C*v.
%% Fortunately C is diagonal.
%M(BcLeftIglob) = M(BcLeftIglob)+half_dt*BcLeftC;
%% Right
%ng = NELY*(NGLL-1)+1;
%BcRightIglob = zeros(ng,1);
%BcRightC = zeros(ng,1);
%for ey=1:NELY,
%  ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
%  e=(ey-1)*NELX+NELX;
%  BcRightIglob(ip) = iglob(NGLL,1:NGLL,e);
%  BcRightC(ip) = BcRightC(ip) + dy_deta*wgll*impedance ;
%end
%M(BcRightIglob) = M(BcRightIglob)+half_dt*BcRightC;
% Top
ng = NELX*(NGLL-1)+1;
BcTopIglob = zeros(ng,1);
BcTopC = zeros(ng,1);
for ex=1:NELX,
  ip = (NGLL-1)*(ex-1)+[1:NGLL] ;
  e=(NELY-1)*NELX+ex;
  BcTopIglob(ip) = iglob(1:NGLL,NGLL,e);
  BcTopC(ip) = BcTopC(ip) + dx_dxi*wgll*impedance ;
end
M(BcTopIglob) = M(BcTopIglob)+half_dt*BcTopC;


%-- DYNAMIC FAULT at bottom boundary
FaultNglob = NELX*(NGLL-1)+1;
FaultIglob = zeros(FaultNglob, 1); 
FaultB = zeros(FaultNglob, 1); 
for ex=1:NELX,
  ip = (NGLL-1)*(ex-1)+[1:NGLL];
  e = ex;
  FaultIglob(ip) = iglob(1:NGLL,1,e);
  FaultB(ip) = FaultB(ip) + dx_dxi*wgll;
end
FaultZ = M(FaultIglob)./FaultB /dt;
FaultX = x(FaultIglob);
FaultV = zeros(FaultNglob,NT);
FaultSigma = repmat(10.,FaultNglob,1);
FaultInitStress = repmat(6.1,FaultNglob,1);
%FaultInitStress( abs(FaultX)<=2 ) = 0.601;
FaultState = zeros(FaultNglob,1);
FaultFriction.MUs = repmat(0.6,FaultNglob,1);
FaultFriction.MUd = repmat(0.5,FaultNglob,1);
%FaultFriction.MUs(abs(FaultX)>10) = 100; % barrier
%FaultFriction.MUd(abs(FaultX)>10) = 100; % barrier
FaultFriction.Dc  = 1;
FaultFriction.W = (FaultFriction.MUs-FaultFriction.MUd)./FaultFriction.Dc;
FaultStrength = FaultSigma.*friction(FaultState,FaultFriction) - FaultInitStress; % strength excess

impedance=1/2;
%v(FaultIglob) = 0.5* max(-FaultStrength,0)/impedance; % initial velocity

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTyseis = [0:10]';		% y coord of receivers
OUTnseis = length(OUTyseis);	% total number of receivers
OUTxseis = zeros(OUTnseis,1);	% x coord of receivers
%********

[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 40;
OUTit = 0;
OUTindx = Init2dSnapshot(iglob);



%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with dissipation
% alpha=1/2, gamma=1
%

for it=1:NT,

 % predictor
  d_mid = d;
  d = d + dt*v +timco1*a; 
  d_mid = 0.5*(d_mid+d); % mid-step displacement:

 % internal forces at mid-step -K*d(t+1/2) 
 % stored in global array 'a'
  a(:) = 0; 
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d_mid(ig);
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
%  a(BcLeftIglob) = a(BcLeftIglob) - BcLeftC .* v(BcLeftIglob);
%  a(BcRightIglob) = a(BcRightIglob) - BcRightC .* v(BcRightIglob);
  a(BcTopIglob) = a(BcTopIglob) - BcTopC .* v(BcTopIglob) ;

 % fault boundary condition: slip weakening
  FaultVFree = v(FaultIglob) + dt*a(FaultIglob)./M(FaultIglob);
  TauStick = FaultZ .*FaultVFree;
  Tau = min(TauStick,FaultStrength); 
  a(FaultIglob) = a(FaultIglob) - FaultB .*Tau;

 % solve for a_new:
  a = a ./M ;

 % update
  v = v + dt*a;
  d = d + timco4*a;

  FaultState = max(2*d(FaultIglob),FaultState);
  FaultStrength = FaultSigma.*friction(FaultState,FaultFriction)-FaultInitStress;
  FaultV(:,it) = 2*v(FaultIglob);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTyseis,time,OUTv);

%    figure(2)
%    Plot2dSnapshot(x,y,v,OUTindx,[0 0.1]);
%    hold on
%    plot(OUTxseis,OUTyseis,'^')
%    hold off

    drawnow

  end

end % ... of time loop

figure(3)
sdt=10;
sdx=1;
surf(time(1:sdt:end),FaultX(1:sdx:end),FaultV(1:sdx:end,1:sdt:end))
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')
