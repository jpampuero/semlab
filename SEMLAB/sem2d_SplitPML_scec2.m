% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% with Perfectly Matched Layers
% and zero initial conditions
% in a structured undeformed grid.
%
% Version 3: domain = rectangular (+symmetry => only right half of the domain)
%            medium = general (heterogeneous)
%            boundaries = 1 fault + 2 PML + 1 free 
%	     time scheme = leapfrog
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
LX=12.5e3;
LY=7.5e3;
%NELX = 50; NELY = 25; P = 8; % polynomial degree
NELX = 25; NELY = 15; P = 4; % polynomial degree

ABC_B = 0;
ABC_R = 2;
ABC_T = 2;
ABC_L = 0;

%********

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x);

RHO = 2670.;
VS  = 3464.;
MU = RHO*VS^2;

ETA = 0.1;
PML_A = 10;
PML_N = 2;

%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
W = wgll * wgll' ;

M     = zeros(nglob,1);		% global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

%**** Set here the parameters of the time solver : ****
%NT = 2200; TT=0; % number of timesteps
NT=0; TT = 12; % total time in seconds
CFL   = 0.6; 			% stability number = CFL_1D / sqrt(2)
%********

dt    = Inf;  			% timestep (set later)

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi ;
coefint2 = jac/dy_deta ;

% FOR EACH ELEMENT ...
for ey=1:NELY, 
for ex=1:NELX, 

  e = (ey-1)*NELX+ex;
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
  rho(:,:) = RHO;
  mu(:,:)  = MU;
%********

 % Diagonal mass matrix
  M(ig) = M(ig) + W .*rho *jac;

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
if ETA, dt=dt/sqrt(1+2*ETA); end
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;
if NT==0, NT=ceil(TT/dt); end

%-- Initialize kinematic fields, stored in global arrays
%sx = zeros(NGLL,NGLL,NEL);
%sy = zeros(NGLL,NGLL,NEL);
f = zeros(nglob,1);
v = zeros(nglob,1);
d = zeros(nglob,1);

time = (1:NT)'*dt;


%-- Absorbing boundaries (first order): 
% The mass matrix needs to be modified at the boundary
% for the IMPLICIT treatment of the term C*v.
% Fortunately C is diagonal.
impedance = RHO*VS;

if ABC_L ==1 % Left
  [BcLC,iBcL] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'L');
  BcLC = impedance*BcLC;
  M(iBcL)  = M(iBcL)  +half_dt*BcLC;
end

if ABC_R ==1 % Right
  [BcRC,iBcR] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'R');
  BcRC = impedance*BcRC;
  M(iBcR) = M(iBcR) +half_dt*BcRC;
end

if ABC_T ==1  % Top
  [BcTC,iBcT] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi,'T');
  BcTC = impedance*BcTC;
  M(iBcT)   = M(iBcT)   +half_dt*BcTC;
end

%--- PML

anyPML = any([ABC_T ABC_R] ==2);
if anyPML

  ePML = [];
  NEL_PML=0;
  if ABC_R>1, ePML= [NELX:NELX:NEL]; end
  if ABC_T>1, ePML= [ ePML [NEL-NELX+1:NEL] ]; end
  ePML = unique(ePML); % +ordered
  NEL_PML = length(ePML);

 %iPML(k) = global node index of the k-th PML node
 %iglob_PML(i,j,ePML) = PML node index of the (i,j,ePML) node

  iglob_PML = iglob(:,:,ePML);
  [iPML,dum,iglob_PML] = unique(iglob_PML(:));
  iglob_PML = reshape(iglob_PML,[NGLL NGLL NEL_PML]);
  nPML = length(iPML);

  axPML=zeros(nPML,1);
  ayPML=zeros(nPML,1);
  xp = x(iPML);
  yp = y(iPML);
  lx = LX-dxe;
  ly = LY-dye;
  if ABC_R, axPML = half_dt*PML_A*VS/dxe * ((xp-lx)/dxe).^PML_N .*(xp>=lx); end
  if ABC_T, ayPML = half_dt*PML_A*VS/dye * ((yp-ly)/dye).^PML_N .*(yp>=ly); end
  clear xp yp

  f_x = zeros(nPML,1);
  f_y = zeros(nPML,1);
  v_x = zeros(nPML,1);
  v_y = zeros(nPML,1);
  s_x = zeros(NGLL,NGLL,NEL_PML);
  s_y = zeros(NGLL,NGLL,NEL_PML);

  % GLL to PML-node connectivity for elements near PML :

  ex1 = 1 + (ABC_L>1); 
  ex2 = NELX - (ABC_R>1);
  ey1 = 1 + (ABC_B>1);
  ey2 = NELY - (ABC_T>1);
  nearPMLconn = zeros(4,NEL); % 1234 = BRTL

  [exPML,eyPML] = ind2sub([NELX NELY], ePML); 
  % right
  ep = find( exPML==NELX & eyPML>=ey1 & eyPML<=ey2 ); 
  nearPMLconn(2,ePML(ep)-1) = ep;
  % top
  ep = find( eyPML==NELY & exPML>=ex1 & exPML<=ex2 ); 
  nearPMLconn(3,ePML(ep)-NELX) = ep;
  
  eNearPML = find( any(nearPMLconn) );
  nearPMLconn = nearPMLconn(:,eNearPML);

  kB = [1:NGLL]';
  kR = [NGLL:NGLL:NGLL*NGLL]';
  kT = kB + NGLL*(NGLL-1);
  kL = kR - NGLL+1;

  nepn = length(eNearPML);
  for epn=1:nepn,
    k=[];
    igp=[];
    if nearPMLconn(2,epn), 
      k=kR; 
      ig=iglob_PML(:,:,nearPMLconn(2,epn)); 
      igp=ig(kL); 
    end
    if nearPMLconn(3,epn), 
      k=[k ; kT]; 
      ig=iglob_PML(:,:,nearPMLconn(3,epn)); 
      igp=[igp ; ig(kB)]; 
    end
    [nearPML(epn).kgll, k] = unique(k);
    nearPML(epn).igPML = igp(k); 
  end

else
  NEL_PML=0;
  ePML =0;
  eNearPML=0;

end


%-- DYNAMIC FAULT at bottom boundary

[FltB,iFlt] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi,'B');
FltN = length(iFlt);
%FltZ = M(iFlt)./FltB /half_dt;  % if friction enforced at n+1/2 ==> dt/2 factor
FltZ = M(iFlt)./FltB /dt;
FltX = x(iFlt);
FltV = zeros(FltN,NT);
FltD = zeros(FltN,NT);
% background stress
FltNormalStress = 120e6;
FltInitStress = repmat(70e6,FltN,1);
% nucleation 
isel = find(abs(FltX)<=1.5e3);
FltInitStress(isel) = 81.6e6;
% friction
FltState = zeros(FltN,1);
FltFriction.MUs = repmat(0.677,FltN,1);
FltFriction.MUd = repmat(0.525,FltN,1);
FltFriction.Dc  = 0.4;
% barrier
L_BARRIER = 15e3/2;
isel = find(abs(FltX)>L_BARRIER);
FltFriction.MUs(isel) = 1e4;
FltFriction.MUd(isel) = 1e4;
FltFriction.W = (FltFriction.MUs-FltFriction.MUd)./FltFriction.Dc;
FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                - FltInitStress; % strength excess

if ETA,  % Kelvin-Voigt viscosity
  NEL_ETA = min( NELX, ceil(L_BARRIER/dxe)+2 );
  x1 = 0.5*(1+xgll');
  eta_taper = exp(-pi*x1.^2); 
  eta = ETA*dt *repmat(eta_taper, NGLL,1 );
else
  NEL_ETA = 0;
end

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = [0:0.5e3:10e3]';		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(2e3,OUTnseis,1);	% y coord of receivers
%********
% receivers are relocated to the nearest node
% OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = floor(0.5/dt);
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);



%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Leapfrog:
%   d[n+1/2] = d[n-1/2] +dt*v[n]
%   v[n+1] = v[n] + dt*( -K*d[n+1/2] +F[n+1/2] )
%
for it=1:NT,

 % update displacement and slip, at mid-step
  d = d +dt*v;
  FltD(:,it) = 2*d(iFlt);
  
 % internal forces -K*d(t+1/2) 
 % stored in global array 'f'
  f(:) = 0;
  f_y(:) = 0;
  f_x(:) = 0;
  ep =1; % PML element counter
  eep = ePML(ep); % next element that is PML
  epn = 1; % near-PML counter
  eepn = eNearPML(epn); % next element that is near PML

  for e=1:NEL,

    isPML = e==eep;
    isNearPML = e==eepn;
    isETA = e<=NEL_ETA;

   %switch to local (element) representation
    ig = iglob(:,:,e);
    if isETA
      local = d(ig) +eta.*v(ig); % Kelvin-Voigt viscosity
    elseif isPML
      local = v(ig);
    else
      local = d(ig);
    end

   %gradients wrt local variables (xi,eta)
   %stress = 2*mu* strain = mu*grad_displ
    sx = MU* Ht*local/dx_dxi ;	
    sy = MU* local*H /dy_deta ;

    if isPML
      igPML = iglob_PML(:,:,ep);
      ax = axPML(igPML);
      sx = ( dt*sx +(1-ax).*s_x(:,:,ep) ) ./(1+ax);
      ay = ayPML(igPML);
      sy = ( dt*sy +(1-ay).*s_y(:,:,ep) ) ./(1+ay);
    end

   %element contribution to internal forces
   %local = coefint1*H*( W.*sx ) + coefint2*( W.*sy )*Ht ;
    d_xi = W.*sx;
    d_xi = coefint1* H * d_xi;
    d_eta = W.*sy;
    d_eta = coefint2* d_eta *Ht;

   %assemble into global vector
    if isPML
      f_x(igPML) = f_x(igPML) - d_xi;
      f_y(igPML) = f_y(igPML) - d_eta;
      s_x(:,:,ep) = sx;
      s_y(:,:,ep) = sy;

      ep = ep+isPML;
      if ep<=NEL_PML, eep = ePML(ep); else eep = 0; end

    else
      local = d_xi  + d_eta ;
      f(ig) = f(ig) -local;

      if isNearPML
        igPML = nearPML(epn).igPML; 
	k = nearPML(epn).kgll;      
        f_x(igPML) = f_x(igPML) - d_xi(k);
        f_y(igPML) = f_y(igPML) - d_eta(k);

        epn = epn+1;
	if epn<=nepn, eepn = eNearPML(epn); else eepn = 0; end
      end

    end

  end 


 % absorbing boundaries:
  if ABC_L==1, f(iBcL) = f(iBcL) - BcLC.*v(iBcL); end
  if ABC_R==1, f(iBcR) = f(iBcR) - BcRC.*v(iBcR); end
  if ABC_T==1, f(iBcT) = f(iBcT) - BcTC.*v(iBcT); end
  if anyPML
    v_x = ( dt*f_x./M(iPML) +(1-axPML).*v_x ) ./(1+axPML);
    v_y = ( dt*f_y./M(iPML) +(1-ayPML).*v_y ) ./(1+ayPML);
  end

 % fault boundary condition: slip weakening
  FltState = max(FltD(:,it),FltState);
  FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                  -FltInitStress;
  %NOTE: half_dt* to enforce friction at mid-step, see also FltZ
  %FltVFree = v(iFlt) + half_dt*f(iFlt)./M(iFlt); 
  FltVFree = v(iFlt) + dt*f(iFlt)./M(iFlt); 
  TauStick = FltZ .*FltVFree;
 % TauStick = f(iFlt)./FltB;
  Tau = min(TauStick,FltStrength); 
  f(iFlt) = f(iFlt) - FltB .*Tau;

 % update velocity
  v = v + dt*f./M;
  if anyPML, v(iPML)= v_x + v_y; end

 % update slip rate
  FltV(:,it) = 2*v(iFlt);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[0 2]); % [0 2]
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off

    drawnow

  end

  test(it) = max(abs(v));
end % ... of time loop

figure(3)
sdt=1;
sdx=4;
PlotSeisTrace(FltX(1:sdx:end),time(1:sdt:end),FltV(1:sdx:end,1:sdt:end));
%surf(time(1:sdt:end),FltX(1:sdx:end),FltV(1:sdx:end,1:sdt:end))
%xlabel('Time')
%ylabel('Position along fault')
%zlabel('Slip rate')
