% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% strain weakening (plasticity) in the bulk,
% paraxial absorbing boundary conditions
% and zero initial conditions
% in a structured undeformed grid.
%
% Version 2d_friction_plastic: 
%	domain = rectangular
%	medium = homogeneous + Coulomb plasticity
%	boundaries = 3 paraxial + 1 slip weakening fault
%	time scheme = central difference
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

% to convert movies to .mov:
MAKEMOVIE = 0;
if MAKEMOVIE & isempty(strfind(path,'makeqtmovie')) 
  addpath ~/local/makeqtmovie/ ; 
end

LX=30; %60;
LY=5; %10;
NELX = 120;
NELY = 20;
P = 4; % polynomial degree
CFL   = 0.45; 	% stability number = CFL_1D / sqrt(2) * EGaOs
		% where EGaOs is in [0.7:1] (=1 if no dissipation)
NT = 2000; %4000; % number of timesteps

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;

NGLL = P+1; % number of GLL nodes per element

[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x);

% Physical properties of the medium
rho = 1;
mu  = 1;
vs = sqrt(mu./rho); 

% initial stress
TauXZInit = 0;
TauYZInit = 1;
% yield for strain weakening plasticity
ys = 1.5;
yd = 0.5;
eplc = 0.2*(ys-yd)/(2*mu)/0.44;
wrpl = (ys-yd)/eplc;
% Maxwell visco-plastic time
tvisc = 0.3*dxe/vs;


%------------------------------------------
% STEP 2: INITIALIZATION

[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
W = wgll * wgll' ;

M     = zeros(nglob,1);		% global mass matrix, diagonal

% plastic parameters
yield = repmat(ys,[NGLL NGLL NEL]);
eplx = zeros(NGLL,NGLL,NEL);
eply = zeros(NGLL,NGLL,NEL);
epl = zeros(NGLL,NGLL,NEL);

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi ;
coefint2 = jac/dy_deta ;

for e=1:NEL, 
  ig = iglob(:,:,e);
  M(ig) = M(ig) + W .*rho *jac;
end

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
dx = diff(xgll)*0.5*min(dxe,dye);
dt = min( dx/vs );
dt = CFL*dt;
half_dt = 0.5*dt;

viscoef = 1+tvisc/dt;


%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;


%-- Absorbing boundaries (first order): 
% M <-- M + 0.5*dt*C
impedance = sqrt(rho*mu);
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
% Right
ng = NELY*(NGLL-1)+1;
BcRightIglob = zeros(ng,1);
BcRightC = zeros(ng,1);
for ey=1:NELY,
  ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
  e=(ey-1)*NELX+NELX;
  BcRightIglob(ip) = iglob(NGLL,1:NGLL,e);
  BcRightC(ip) = BcRightC(ip) + dy_deta*wgll*impedance ;
end
M(BcRightIglob) = M(BcRightIglob)+half_dt*BcRightC;
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
FaultZ = M(FaultIglob)./(half_dt*FaultB);
FaultX = x(FaultIglob);
FaultV = zeros(FaultNglob,NT);
FaultRuptureFront = zeros(FaultNglob,1);
%**** Set here friction properties and initial conditions for the fault : ****
LF = 40; % final half-size
FaultState = zeros(FaultNglob,1);
FaultFriction.MUs = repmat(1.2,FaultNglob,1);
FaultFriction.MUd = repmat(0.8,FaultNglob,1);
FaultFriction.MUs(abs(FaultX)>LF) = 100; % barrier
FaultFriction.MUd(abs(FaultX)>LF) = 100; % barrier
FaultFriction.Dc  = 1;
FaultFriction.W = (FaultFriction.MUs-FaultFriction.MUd)./FaultFriction.Dc;
Lc = 2*0.57888694 *mu/FaultFriction.W(1) ; % nucleation size
Sm = FaultFriction.W(1)/(mu/2/vs);
FaultInitStress = zeros(FaultNglob,1);
ix = find(abs(FaultX)<Lc/2);
FaultInitStress(ix) = 1.2+1e-4;
ix = find(abs(FaultX)>=Lc/2);
FaultInitStress(ix) = 1. +0.2* exp(-(abs(FaultX(ix))-Lc/2).^2/(4*Lc)^2 ) ;
%FaultInitStress = 6 +1e-4* exp(-(abs(FaultX)-Lc/2).^2/(0.4*Lc)^2 ) ;
%********

%-- initialize data for output seismograms
OUTxseis = [0:0.05:1]'*LX;
OUTnseis = length(OUTxseis);	% total number of receivers
OUTyseis = repmat(LY/2, OUTnseis,1);
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 200; %400;
OUTindx = Plot2dSnapshot(iglob);
OUTit = 0;

if MAKEMOVIE
  MakeQTMovie('start', 'plastic.mov')
  MakeQTMovie('framerate', 12)
  MakeQTMovie('quality', .8)
end
  

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit central difference
%
% LINEAR STRAIN WEAKENING PLASTICITY (Andrews 1976)
%   yield(epl) = max( ys-wrpl*epl , yd )
% return map algorithm (Wilkins 1964, see Zienkiewicz and Taylor, Vol.2)
% with (linear) PERZYNA VISCOPLASTIC REGULARIZATION
%   depl/dt = 1/tvisc *max( |stress| - yield, 0 )/(2*mu)

for it=1:NT,

  v = v+ half_dt*a;	% = v_tmp, partial velocity update
  d = d + dt*v; 	% explicit displacement update 
  a(:) = 0; 		% forces are temporarily stored in global array 'a'

 % friction update
  FaultState = max(2*d(FaultIglob),FaultState);
  FaultStrength = friction(FaultState,FaultFriction)-FaultInitStress;

  for e=1:NEL,

 % STEP 1: TRIAL STRESS (elastic update)
 % internal forces at new step -K*d(t+1) 

   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);

   %gradients wrt local variables (xi,eta) = Ht*local and local*H
   %stress = 2*mu*( strain - plastic_strain )
    sx = mu * ( Ht*local/dx_dxi -2*eplx(:,:,e) );	
    sy = mu * ( local*H/dy_deta -2*eply(:,:,e) );
 
   %test yield
    sxa = TauXZInit +sx;
    sya = TauYZInit +sy;
    sabs = sqrt( sxa.^2 + sya.^2 );
    phi = sabs-yield(:,:,e);

 % STEP 2: RETURN MAPPING (plastic solver)
    if any(phi(:)>0)
      % solve for plastic |strain| increment, depl in
      %   depl/dt = 1/tvisc*( |stress_trial| -2*mu*depl - yield(epl+depl) )/(2*mu)
      % with linear strain weakening:
      %   yield(epl) = max( ys-wrpl*epl , yd )
      %   yield(epl+depl) = max( yield(epl)-wrpl*depl, yd )
      depl = min( phi/(viscoef*2*mu - wrpl), (sabs - yd)/(viscoef*2*mu) );
      depl = max(depl,0);
      % update yield
      yield(:,:,e) = max( yield(:,:,e)-wrpl*depl, yd );
      % update stress
      nx = sxa./sabs;
      ny = sya./sabs;
      ds = -2*mu*depl;
      sx = sx + ds.*nx;
      sy = sy + ds.*ny;
      % update plastic strain
      epl(:,:,e) = epl(:,:,e) + depl;
      eplx(:,:,e) = eplx(:,:,e) + depl.*nx;
      eply(:,:,e) = eply(:,:,e) + depl.*ny;
    end

   %element contribution to internal forces
   %local = coefint1*H*( W.*sx ) + coefint2*( W.*sy )*Ht ;
    d_xi = W.*sx;
    d_xi = H * d_xi;
    d_eta = W.*sy;
    d_eta = d_eta *Ht;
    local = coefint1* d_xi  + coefint2* d_eta ;
   %assemble into global vector
    a(ig) = a(ig) -local;

  end 

 % boundary conditions
 % implicit absorbing boundaries: 
 % here we use v_tmp, the rest is done through the correction of M
%  a(BcLeftIglob) = a(BcLeftIglob) - BcLeftC .* v(BcLeftIglob);
  a(BcRightIglob) = a(BcRightIglob) - BcRightC .* v(BcRightIglob);
  a(BcTopIglob) = a(BcTopIglob) - BcTopC .* v(BcTopIglob) ;

 % fault boundary condition: slip weakening
  FaultVFree = v(FaultIglob) + half_dt*a(FaultIglob)./M(FaultIglob);
  TauStick = FaultZ .*FaultVFree;
 % TauStick = a(FaultIglob)./FaultB;
  Tau = min(TauStick,FaultStrength); 
  a(FaultIglob) = a(FaultIglob) - FaultB .*Tau;

 % solve for a_new
  a = a ./M ;

 % velocity update: v_new = v_old + dt/2*( a_old+a_new )
  v = v + half_dt*a;

  FaultV(:,it) = 2*v(FaultIglob);
  FaultRuptureFront( FaultRuptureFront==0 & FaultV(:,it)>1e-3 ) = it*dt;

%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);

  if mod(it,OUTdt) == 0

    OUTit = OUTit+1; 

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);
    ylabel('X')

    figure(2)
    clf
    subplot(211)
    Plot2dSnapshot(x,y,v,OUTindx); %,[-1 1]);
    title('Velocity')
    xlabel('')
    hold on; plot(OUTxseis,OUTyseis,'^'); hold off
    %rect = get(gcf,'Position');
    %rect(1:2) = [0 0];
    %OUTmovie1(:,OUTit)=getframe(gcf,rect);

%    figure(3)
    subplot(212)
    eplm(iglob) = epl;
    Plot2dSnapshot(x,y,eplm/eplc,OUTindx); % ,[0 10]);
    title('Plastic strain / \epsilon_c')
    %rect = get(gcf,'Position');
    %rect(1:2) = [0 0];
    %OUTmovie2(:,OUTit)=getframe(gcf,rect);
    drawnow

    if MAKEMOVIE, MakeQTMovie('addfigure'); end

  end

end % ... of time loop

if MAKEMOVIE, MakeQTMovie('finish'); end
