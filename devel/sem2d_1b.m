% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source,
% in a structured undeformed grid.
%
% Version 1b:	domain = rectangular
%            	medium = general (heterogeneous)
%		no global arrays
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

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

LX=10;
LY=30;
NELX = 10;
NELY = 30;
dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;



%------------------------------------------
% STEP 2: INITIALIZATION

P = 8; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 1200; % number of timesteps

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

W     = zeros(NGLL,NGLL,NELX,NELY);	% for internal forces
x     = zeros(NGLL,NGLL,NELX,NELY);	% coordinates of GLL nodes
y     = zeros(NGLL,NGLL,NELX,NELY);	
M     = zeros(NGLL,NGLL,NELX,NELY);	% global mass matrix, diagonal
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

 % Coordinates of the computational (GLL) nodes
  x(:,:,e) = repmat( dxe*(ex-0.5+0.5*xgll) , 1,NGLL);
  y(:,:,e) = repmat( dye*(ey-0.5+0.5*xgll'), NGLL,1);

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  rho(:,:) = 1;
  if ex>NELX/2-2 & ex<NELX/2+3
    mu(:,:)  = 0.25;
  else
    mu(:,:)  = 1;
  end

 % The diagonal mass matrix is stored in a local array. 
 % It will be assembled later
 % (nodes at the boundary between two elements get 
 % contributions from both).
  M(:,:,e) = M(:,:,e) + wgll2 .*rho *jac;

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

%-- Assemble mass matrix
% Along X
M(1,:,2:NELX,:) = M(1,:,2:NELX,:) + M(NGLL,:,1:NELX-1,:);
M(NGLL,:,1:NELX-1,:) = M(1,:,2:NELX,:);
% Along Y
M(:,1,:,2:NELY) = M(:,1,:,2:NELY) + M(:,NGLL,:,1:NELY-1);
M(:,NGLL,:,1:NELY-1) = M(:,1,:,2:NELY);


%-- Initialize kinematic fields, stored in local arrays
d = zeros(NGLL,NGLL,NELX,NELY);
v = zeros(NGLL,NGLL,NELX,NELY);
a = zeros(NGLL,NGLL,NELX,NELY);

time = (1:NT)'*dt;

%-- SOURCE TERM: point force, time function = Ricker wavelet
% located in element (ex,ey), local node (igll,jgll)
% NOTE: in this example sources at element vertex (connectivity = 4)
ex = [NELX/2;NELX/2+1;NELX/2;  NELX/2+1] ;
ey = [NELY/2;NELY/2;  NELY/2+1;NELY/2+1] ;
igll = [NGLL;1   ;NGLL;1];
jgll = [NGLL;NGLL;1   ;1];
Fig = sub2ind([NGLL NGLL NELX NELY],igll,jgll,ex,ey);
Ff0 = 0.3; 	% fundamental frequency
Ft0 = 1.5/Ff0; 	% delay
% source time function (at mid-steps)
Ft = ricker( time-0.5*dt, Ff0,Ft0);

%-- initialize data for output seismograms
% record at these nodes:
ey = repmat( (NELY/2-7:NELY/2+8), 2,1); ey=ey(:);
OUTnseis = length(ey);
ex = repmat( NELX/2-2, OUTnseis,1);
igll = repmat(P/2+1, OUTnseis,1);
jgll = repmat([1;P/2+1], OUTnseis/2,1);
OUTiglob = sub2ind([NGLL NGLL NELX NELY],igll,jgll,ex,ey);
OUTxseis = x(OUTiglob); 
OUTyseis = y(OUTiglob); 
OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTdt = 50;
OUTit = 0;
OUTindx = Init2dSnapshot_b(NGLL);


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

 % internal forces at mid-step -K*d(t+1/2) 
 % stored in array 'a'
  for e=1:NEL, 
   % NOTE: assuming a(:,:,ex,ey) morphs into a(:,:,e)
   %gradients wrt local variables (xi,eta)
    local = d(:,:,e);
    d_xi  = Ht*local;	
    d_eta = local*H;
   %element contribution to internal forces
   %a(:,:,e) = coefint1*H*( W(:,:,e).*d_xi )+ coefint2*( W(:,:,e).*d_eta )*Ht;
    wloc = W(:,:,e);
    d_xi = wloc.*d_xi;
    d_xi = H * d_xi;
    d_eta = wloc.*d_eta;
    d_eta = d_eta *Ht;
    a(:,:,e) = coefint1* d_xi  + coefint2* d_eta ;
  end 

 % assemble forces
 % Along X
  a(1,:,2:NELX,:) = a(1,:,2:NELX,:) + a(NGLL,:,1:NELX-1,:);
  a(NGLL,:,1:NELX-1,:) = a(1,:,2:NELX,:);
 % Along Y
  a(:,1,:,2:NELY) = a(:,1,:,2:NELY) + a(:,NGLL,:,1:NELY-1);
  a(:,NGLL,:,1:NELY-1) = a(:,1,:,2:NELY);

 % add external forces
 % NOTE: in Matlab, a(Fig) = a(Figll,Fjgll,Fex,Fey)
  a(Fig) = a(Fig) - Ft(it);

 % acceleration: a = (-K*d +F)/M
  a = - a ./M ;

 % update
 % v_new = v_old + dt*a_new;
 % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
 %       = d_mid + 0.5*dt*v_new
  v = v + dt*a;
  d = d + half_dt*v;


%------------------------------------------
% STEP 4: OUTPUT

  OUTv(:,it) = v(OUTiglob);
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(3) % seismograms
    PlotSeisTrace(OUTyseis,time,OUTv);

    figure(2)
    Plot2dSnapshot_b(x,y,v,OUTindx,0.5);
    hold on
    plot(OUTxseis,OUTyseis,'^',x(Fig),y(Fig),'*')
    hold off

    drawnow

  end

end % ... of time loop
