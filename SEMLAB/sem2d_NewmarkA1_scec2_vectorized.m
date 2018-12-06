% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% dynamic fault with slip weakening,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid.
%
% You can modify simulation parameters in parts of the code 
% that have comments starting as: "**** Set here ..."
%
% Version 2: domain = rectangular
%            medium = general (heterogeneous)
%            boundaries = 1 fault + 3 paraxial
%	     time scheme = Newmark as in SPECFEM3D: alpha=1,beta=0,gamma=1/2
%
% Aug 7 2007: + effect of finite seismogenic depth LZ (2.5D, crude crustal-plane model)
%             + Kelvin-Voigt viscosity (ETA)
%             + moved friction update to begining of timestep
%             + immediate healing when slip rate = 0
%
% May 29 2018: + optional symmetry with respect to x=0
%            	 (then boundaries = 1 fault + 2 paraxial + 1 free)
%
% July 13 2018: vectorized version, runs faster
%               needs multiprod from 
%               https://fr.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	ampuero@gps.caltech.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

disp('Initializing ...')

%**** Set here the parameters of the square box domain and mesh : ****
LX=50e3;
LY=50e3/3;
LZ=10e3;  % thickness of the seismogenic region (approximately accounted for). Turned off if LZ=inf.
%NELX = 75; NELY = 25; P = 8; % polynomial degree
NELX = 150; NELY = 50; P = 4; % polynomial degree
SYM_X = 1; % If SYM_X=1, enforce symmetry with respect to x=0
           % (the left boundary becomes a free surface)
%********

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
if ~SYM_X, x = x-LX/2; end
nglob = length(x);

RHO = 2670.;
VS  = 3464.;

ETA = 0.1; % Kelvin-Voigt viscosity term = ETA*dt*K*v

%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

isCrustalPlane = isfinite(LZ);

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
M     = zeros(nglob,1);		% global mass matrix, diagonal
if (isCrustalPlane), KZ = zeros(nglob,1); end % global crustal-plane matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

%**** Set here the parameters of the time solver : ****
NT = 1000; %2200; % number of timesteps
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
for e=1:NEL,
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
  rho(:,:) = RHO;
  mu(:,:)  = RHO* VS^2;
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

 % diagonal matrix for crustal plane term, KZ*d = mu*(pi/2/LZ)^2 * displacement
 % pi/2 comes from quarter-of-wavelength proxy
  if (isCrustalPlane),
    KZ(ig) = KZ(ig) + wgll2 .*mu *jac *(pi/2/LZ)^2 ;
  end

end %... of element loop

dt = CFL*dt;
if ETA, dt=dt/sqrt(1+2*ETA); end
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;

% connectivity matrix for vectorized assembly (without loops)
Conn = sparse(iglob(:),[1:NGLL*NGLL*NEL],1);

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

a_elem = zeros(NGLL,NGLL,NEL);

time = (1:NT)'*dt;

%-- Absorbing boundaries (first order): 
% The mass matrix needs to be modified at the boundary
% for the IMPLICIT treatment of the term C*v.
% Fortunately C is diagonal.
impedance = RHO*VS;

if ~SYM_X,
  [BcLC,iBcL] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'L'); % Left
  BcLC = impedance*BcLC;
  M(iBcL) = M(iBcL) +half_dt*BcLC;
end

[BcRC,iBcR] = BoundaryMatrix(wgll,[NELX NELY],iglob,dy_deta,'R'); % Right
BcRC = impedance*BcRC;
M(iBcR) = M(iBcR) +half_dt*BcRC;

[BcTC,iBcT] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi ,'T'); % Top
BcTC = impedance*BcTC;
M(iBcT) = M(iBcT) +half_dt*BcTC;

%-- DYNAMIC FAULT at bottom boundary
[FltB,iFlt,jFlt] = BoundaryMatrix(wgll,[NELX NELY],iglob,dx_dxi,'B');
FltN = length(iFlt);
FltZ = M(iFlt)./FltB /half_dt;
FltX = x(iFlt);
FltV = zeros(FltN,NT+1);
FltD = zeros(FltN,NT+1);
% background
FltNormalStress = 120e6;
FltInitStress = repmat(70e6,FltN,1);
FltState = zeros(FltN,1);
FltFriction.MUs = repmat(0.677,FltN,1);
FltFriction.MUd = repmat(0.525,FltN,1);
FltFriction.Dc  = 0.4;
% barrier
L_BARRIER = 15e3; % hypocentral distance to barrier
isel = find(abs(FltX)>L_BARRIER);
FltFriction.MUs(isel) = 1e4; % barrier
FltFriction.MUd(isel) = 1e4; % barrier
% nucleation
isel = find(abs(FltX)<=1.5e3);
FltInitStress(isel) = 81.6e6;
FltFriction.W = (FltFriction.MUs-FltFriction.MUd)./FltFriction.Dc;
FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                - FltInitStress; % strength excess

if ETA,  % Kelvin-Voigt viscosity
  isKelvinVoigt = zeros(NEL,1);
  if SYM_X,
    e = find( [dxe/2:dxe:(LX-dxe/2)]<L_BARRIER+dxe );
  else
    e = find( abs([-(LX/2+dxe/2):dxe:(LX/2-dxe/2)])<L_BARRIER+dxe );
  end
  isKelvinVoigt(e) = 1;
  x1 = 0.5*(1+xgll');
  eta_taper = exp(-pi*x1.^2); 
  eta_elem = ETA*dt *repmat(eta_taper, NGLL,1 );
  eta = zeros(NGLL,NGLL,NEL);
  for e =1:NEL,
    if isKelvinVoigt(e), eta(:,:,e) = eta_elem; end
  end
end

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = [-16e3:600:16e3]';		% x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(7.5e3,OUTnseis,1);	% y coord of receivers
%********
% receivers are relocated to the nearest node
% OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 50;
OUTit = 0;
OUTindx = Plot2dSnapshot(iglob);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark scheme with
% alpha=1, beta=0, gamma=1/2
%
disp('Starting time loop ...')

for it=1:NT,


 % update
  d = d + dt*v + half_dt_sq*a; 
%  FltState = max(2*d(iFlt),FltState); % no healing
  FltState = ( FltState + max(2*d(iFlt)-FltD(:,it), 0) ) .* (v(iFlt)>0); % with healing

 % prediction 
  v = v + half_dt*a;

 % internal forces -K*d(t+1) 
 % stored in global array 'a'
 %switch to local representation (element-by-element) 
  a_elem = d(iglob) + eta.*v(iglob);

 % vectorized version (without loops) using multiprod from 
 % https://fr.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications--with-array-expansion-enabled
 %gradients wrt local variables (xi,eta)
  d_xi  = multiprod( Ht , a_elem );
  d_eta = multiprod( a_elem, H );
 %element contribution to internal forces
 %local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
  d_xi = W .* d_xi;
  d_xi = multiprod( H , d_xi );
  d_eta = W .* d_eta;
  d_eta = multiprod( d_eta , Ht);
  a_elem = coefint1 * d_xi  + coefint2 * d_eta ;

 %assemble into global vector
 %vectorized (without loop) using a sparse connectivity matrix:
  a = - Conn * a_elem(:);

 % absorbing boundaries:
  if ~SYM_X,  a(iBcL) = a(iBcL) - BcLC .* v(iBcL); end
  a(iBcR) = a(iBcR) - BcRC .* v(iBcR);
  a(iBcT) = a(iBcT) - BcTC .* v(iBcT) ;

 %crustal-plane term
  if (isCrustalPlane), a = a - KZ.*d; end

 % fault boundary condition: slip weakening
  FltStrength = friction(FltState,FltFriction)*FltNormalStress ...
                  -FltInitStress;
%% rupture nucleation through time weakening, like Andrews 76
%  VNUC = 0.4;
%%  if time(it)<10/VNUC
%    %FltStrength = min( FltStrength,...
%    ix = find(abs(FltX)<=10);
%    FltStrength(ix) = min(FltStrength(ix), ...
%              max(0.5,0.55+(abs(FltX(ix))-VNUC*time(it))*0.05/(1.0*dxe))...
%                     - FltInitStress(ix) ) ;
%%  end

  FltVFree = v(iFlt) + half_dt*a(iFlt)./M(iFlt);
  TauStick = FltZ .*FltVFree;
 % TauStick = a(iFlt)./FltB;
  Tau = min(TauStick,FltStrength); 
  a(iFlt) = a(iFlt) - FltB .*Tau;

 % solve for a_new:
  a = a ./M ;

 % correction
  v = v + half_dt*a;

  FltV(:,it+1) = 2*v(iFlt);
  FltD(:,it+1) = 2*d(iFlt);

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[0 2]);
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off

    drawnow

  end

end % ... of time loop

%keyboard
disp('Displaying results ...')

figure(3)
subplot(211)
sdt=10;
sdx=4;
surf(time(1:sdt:end),FltX(1:sdx:end),FltV(1:sdx:end,2:sdt:end))
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')


%-- interp1_sem example --
% Interpolate a fault field (slip rate) on a regular grid

% Note that BoundaryMatrix was called with a third output: 
% "jFlt" is a table of indices [GLL node, fault element] --> fault node index.

xi = linspace(FltX(1),FltX(end),length(FltX));	% a regular interpolation grid
db=interp1_sem(FltX,xi,jFlt); 			% initialize the interpolation
for it=1:ceil(NT/sdt),
  vi(:,it) = interp1_sem(FltV(:,2+(it-1)*sdt), db,jFlt); % do the interpolation
end

subplot(212)
surf(time(1:sdt:end),xi,vi)
xlabel('Time')
ylabel('Position along fault')
zlabel('Slip rate')

figure(4)
plot(FltX,FltV(:,300), xi,interp1_sem(FltV(:,300), db,jFlt));
ylabel('X')
ylabel('V')
legend('SEM','interpolated')
