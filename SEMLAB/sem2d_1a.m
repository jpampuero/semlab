% SEM2D	applies the Spectral Element Method
% to solve the 2D SH wave equation, 
% with stress free boundary conditions,
% zero initial conditions
% and a time dependent force source,
% in a structured undeformed grid.
%
% Version 1a:	domain = rectangular
%            	medium = general (heterogeneous)
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
LX=10; 	% x-size of the box
LY=30;
NELX = 10;
NELY = 30;
P = 8; % polynomial degree
%********

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
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
NT = 1200; % number of timesteps
CFL   = 0.6; 			% stability number
%********

dt = inf; % will be set later

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
% . set physical properties
% . set mass and stiffness matrices
% . set timestep
for ey=1:NELY, 
for ex=1:NELX, 

  e = (ey-1)*NELX+ex;
  ig = iglob(:,:,e);

%**** Set here the physical properties of the heterogeneous medium : ****
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
 % example: a low velocity layer
  rho(:,:) = 1;
  if ex>NELX/2-2 & ex<NELX/2+3
    mu(:,:)  = 0.25;
  else
    mu(:,:)  = 1;
  end
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

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

time = (1:NT)'*dt;

%-- SOURCE TERM: point force, time function = Ricker wavelet
%**** Set here the source location : ****
Fx = 5; Fy = 15;
%********
[Fx,Fy,Fig] = FindNearestNode(Fx,Fy,x,y);
Ff0 = 0.3; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = src_timef( time-0.5*dt,'ricker', Ff0,Ft0);

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTyseis = [7.5:0.5:22.5]';		% y coord of receivers
OUTnseis = length(OUTyseis);		% total number of receivers
OUTxseis = repmat(2.5,OUTnseis,1);	% x coord of receivers
%********

[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
OUTv = zeros(OUTnseis,NT);

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
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;

    figure(1) % seismograms
    PlotSeisTrace(OUTyseis,time,OUTv);

    figure(2)
    Plot2dSnapshot(x,y,v,OUTindx,[-0.5 0.5]);
    hold on
    plot(OUTxseis,OUTyseis,'^',Fx,Fy,'*')
    hold off
    
%    rect = get(gcf,'Position'); rect(1:2) = [0 0]; OUTmovie(:,OUTit)=getframe(gcf,rect);

    drawnow

  end

end % ... of time loop
%disp('To replay the movie: movie(OUTmovie)')
