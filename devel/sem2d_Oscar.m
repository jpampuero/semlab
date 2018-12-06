% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% in a structured undeformed grid,
% with stress free conditions on the vertical boundaries,
% paraxial absorbing condition on the bottom boundary,
% zero initial conditions
% and a normally incident plane wave as source.
%
% Version for Oscar Ishizawa: "city effect" (buildings)
% Version 1a:	domain = rectangular
%            	medium = general (heterogeneous)
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

% Typical frame building properties:
% C ~ 100 m/s
% rho ~ 0.05 *rho_soil
% plastic strain limit ~ 0.005

%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,LX]*[0,LY] is divided 
% into NELX*NELY quadrangular elements.
% The numbering of elements follows this convention:
%
%      ... ... ... ... NELX*NELY
% ^    ... ... ... ... ...
% | NELX+1 ... ... ... 2*NELX
% |     1   2  ... ... NELX  
% --->
%
% Actually, in this case the macro-mesh is so simple 
% that we don't need to build and store its database here.
% The SEM mesh (GLL nodes for each element) is built in the next step.

LX=40;
LY=120;
NELX = 2;
NELY = 6;
dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;



%------------------------------------------
% STEP 2: INITIALIZATION

P = 8; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 15000; % number of timesteps

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

iglob = zeros(NGLL,NGLL,NEL);	% local to global index mapping
W     = zeros(NGLL,NGLL,NEL);	% for internal forces
nglob = (NELX*(NGLL-1)+1)*(NELY*(NGLL-1)+1);	% number of global nodes
x     = zeros(nglob,1);		% coordinates of GLL nodes
y     = zeros(nglob,1);	
M     = zeros(nglob,1);		% global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	% density will not be stored
mu    = zeros(NGLL,NGLL);	% shear modulus will not be stored

CFL   = 0.6; 			% stability number
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

 % The table I = iglob(i,j,e) maps the local numbering of the 
 % computational nodes (i,j,e) to their global numbering I.
 % 'iglob' is used to build global data 
 % from local data (contributions from each element)
  if e==1  % first element: bottom-left
    ig = reshape([1:NGLL*NGLL],NGLL,NGLL);
  else
    if ey==1 	%  elements on bottom boundary
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(2:NGLL,:) = iglob(NGLL,NGLL,e-1) ... 	% the rest
                   + reshape([1:NGLL*(NGLL-1)],NGLL-1,NGLL);
    elseif ex==1 % elements on left boundary 
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(:,2:NGLL) = iglob(NGLL,NGLL,e-1) ... 	% the rest
                      + reshape([1:NGLL*(NGLL-1)],NGLL,NGLL-1);
    else 	% other elements
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(2:NGLL,2:NGLL) = iglob(NGLL,NGLL,e-1) ...
                      + reshape([1:(NGLL-1)*(NGLL-1)],NGLL-1,NGLL-1);
    end
  end
  iglob(:,:,e) = ig;

 % Coordinates of the computational (GLL) nodes
  x(ig) = repmat( dxe*(ex-0.5+0.5*xgll) , 1,NGLL);
  y(ig) = repmat( dye*(ey-0.5+0.5*xgll'), NGLL,1);

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  if ey>2*NELY/3
    rho(:,:) = 1400;
    mu(:,:)  = 1400*65^2; % vs=65
  else
    rho(:,:) = 2000;
    mu(:,:)  = 2000*600^2; % vs=600
  end
  
 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
  M(ig) = M(ig) + wgll2 .*rho *jac;

 % The stiffness matrix K is not assembled at this point
 % We only store its local contributions
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

disp(sprintf('Timestep DeltaT = %f',dt))
disp(sprintf('Total duration  = %f',dt*NT))

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;

%-- SOURCE : plane wave with normal incidence, 
% time function = Ricker wavelet
Ff0 = 1.2; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = ricker( time-0.5*dt, Ff0,Ft0);

%-- initialize data for output seismograms
% record at these global nodes:
ex=(1:NELX);
ey=NELY;
igll=[1:NGLL-1]; %[1 P/2+1];
jgll=NGLL;
e=(ey-1)*NELX+ex;
stuff = iglob(igll,jgll,e);   OUTiglob = stuff(1:end)';
OUTxseis = x(OUTiglob); 
OUTyseis = y(OUTiglob); 
OUTnseis = length(OUTiglob);
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 100;
OUTit = 0;
OUTindx = Init2dSnapshot(iglob);

%-- BUILDINGS are modelled as single-degree-of-freedom oscillators
% tied to some surface nodes
ex=(1:NELX);
ey=NELY;
igll= P/2+1;
jgll=NGLL;
e=(ey-1)*NELX+ex;
stuff = iglob(igll,jgll,e);   
iglobB = stuff(:);
MB = 350*100*8; % 100 m^2 per floor, 8 floors
KB = 16e6; % such that f0 ~ 1.2 Hz
nB = length(iglobB);
dB = zeros(nB,1);
vB = zeros(nB,1);
aB = zeros(nB,1);
OUTvB = zeros(nB,NT);
xB = x(iglobB);


%-- Absorbing condition (first order) at bottom boundary :
impedance = 2000*600; % = rho*vs
ng = NELX*(NGLL-1)+1;
BcBottomIglob = zeros(ng,1);
BcBottomC = zeros(ng,1);
for ex=1:NELX,
  ip = (NGLL-1)*(ex-1)+[1:NGLL] ;
  BcBottomIglob(ip) = iglob(1:NGLL,1,ex);
  BcBottomC(ip) = BcBottomC(ip) + dx_dxi*wgll*impedance ;
end
M(BcBottomIglob) = M(BcBottomIglob)+half_dt*BcBottomC;



%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%

for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 
  dB = dB + half_dt*vB; 

 % internal forces at mid-step -K*d(t+1/2) :
  a(:) = 0; % store -K*d in a global array
  for e=1:NEL,
   %switch to local (element) representation
    ig = iglob(:,:,e);
    local = d(ig);
   %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
   %local = coefint1* H * ( W(:,:,e).*d_xi ) + coefint2* ( W(:,:,e).*d_eta ) *Ht ;
    wloc = W(:,:,e);
    d_xi = wloc.*d_xi;
    d_xi = H * d_xi;
    d_eta = wloc.*d_eta;
    d_eta = d_eta *Ht;
    local = coefint1* d_xi  + coefint2* d_eta ;
   %assemble into global vector
    a(ig) = a(ig) -local;
  end 
  aB = -KB.*( dB-d(iglobB) );
  a(iglobB) = a(iglobB) -aB;

 % absorbing boundaries and incident plane wave Ft:
  a(BcBottomIglob) = a(BcBottomIglob) - BcBottomC .* (v(BcBottomIglob)-Ft(it));

 % acceleration: a = (-K*d +F)/M
  a = a ./M ;
  aB = aB ./MB ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;
  vB = vB + dt*aB; 
  dB = dB + 0.5*dt*vB;


%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  OUTvB(:,it) = vB;
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    clf
    subplot(211)
    PlotSeisTrace(OUTxseis,time,OUTv);
    subplot(212)
    PlotSeisTrace(xB,time,OUTvB);
    ylabel('Velocity on buildings')

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[-1 1]);
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off
    
    rect = get(gcf,'Position');
    rect(1:2) = [0 0];
    OUTmovie(:,OUTit)=getframe(gcf,rect);

    drawnow

  end

end % ... of time loop
%disp('To replay the movie: movie(OUTmovie)')
