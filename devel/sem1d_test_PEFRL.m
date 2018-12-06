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

% THIS VERSION: applied force on the bottom boundary
%		absorbing condition on the top boundary

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=10; NEL = 10;
Fx = 0;		% source location
OUTx  = (0:L)'; % receiver locations
ABSO_TOP = 0; ABSO_BOTTOM =0; 	% absorbing boundaries

% the following variables can be set prior to calling the script:
if ~exist('P','var'), P = 8; end % polynomial degree
if ~exist('Ff0','var'), Ff0 = 0.5; end 	% fundamental frequency of the source
if ~exist('CFL','var'), CFL = 1.25; end	% stability number <=1.25
if ~exist('TMAX','var'), TMAX = 20*2*L; end	% total simulation time (2*L= round trip travel time)


%------------------------------------------
% STEP 2: INITIALIZATION

X = [0:NEL]'*L/NEL;
NGLL = P+1; % number of GLL nodes per element
[iglob,coor,nglob] = mesh1d(X,NGLL);	

% Physical properties of the medium
rho = 1;
mu = 1;
[M,K] = BuildMK_1d(coor,iglob,rho,mu);
K = K(:,:,1);

dt0 = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
%*** make travel time over 2*L a multiple of dt ***
dt = 2*L/ceil(2*L/dt0);
disp(sprintf('Changed CFL from %f to %f',CFL,CFL*dt/dt0))
CFL = CFL*dt/dt0;
NT = ceil(TMAX/dt); % number of timesteps

half_dt = 0.5*dt;

if ABSO_TOP
  BcTopC = sqrt(rho*mu);
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho*mu);
  M(1) = M(1)+half_dt*BcBottomC;
end

d = zeros(nglob,1);
v = zeros(nglob,1);
f = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
[Fdist,Fix] = min( abs(coor-Fx) );
% source time function 
Ft0 = 1.5/Ff0; % delay
t = (1:NT)'*dt;
Ft = src_timef(t,'ricker', Ff0,Ft0);
%Ft = src_timef(t,'gabor', Ff0,[], 3 );
%Ft = src_timef(t,'gaussian', Ff0,Ft0 );

% output arrays
OUTnx = length(OUTx);
OUTix = zeros(OUTnx,1);
% relocate to nearest GLL node
OUTdist = zeros(OUTnx,1);
for i=1:OUTnx,
  [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) );
end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,NT);
OUTv = zeros(OUTnx,NT);

%v = ricker( coor, Ff0,Ft0);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% PEFRL time scheme (Omelyan, Mryglod and Folk, 2002)
%

xi = 0.1786178958448091E+00;
lambda = -0.2123418310626054E+00;
chi = -0.6626458266981849E-01;

coefd = [xi; chi; 1-2*(chi+xi); chi; xi] *dt;
coeff = cumsum(coefd(1:4));
coefv = [(1-2*lambda)/2; lambda; lambda; (1-2*lambda)/2] *dt;

for it=1:NT,
 Ftl = src_timef((it-1)*dt+coeff,'ricker',Ff0,Ft0);
 for nstage=1:4, 

  d = d + coefd(nstage)*v;

 % internal forces
  f(:) = 0;
  for e=1:NEL,
    ix = iglob(:,e);
    f(ix) = f(ix) - K*d(ix) ;
  end 
 % add external forces
  %f(Fix) = f(Fix) + Ft(it);
  f(Fix) = f(Fix) + Ftl(nstage);

  v = v +coefv(nstage)*f./M ;

 end % ... of the 4 stages
  d = d + coefd(5)*v;

%------------------------------------------
% STEP 4: OUTPUT
  
  OUTd(:,it) = d(OUTix);
  OUTv(:,it) = v(OUTix);

end % of time loop
