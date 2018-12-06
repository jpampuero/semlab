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
F_IS_WAVE = 0;	% 1= incident wave from bottom
OUTx  = (0:L)'; % receiver locations
ABSO_TOP = 0; ABSO_BOTTOM =0; 	% absorbing boundaries
NEL_VISC = 0; tvisc = 0.02;	% viscous layer

% the following variables can be set prior to calling the script:
if ~exist('CFL','var'), CFL = 0.85 ; end %/2.7; % stability number <=0.85
if ~exist('TMAX','var'), TMAX = 20*2*L;	end % total simulation time (2*L= round trip travel time)
if ~exist('P','var'), P = 8; end % polynomial degree
if ~exist('Ff0','var'), Ff0 = 0.5; end 	% fundamental frequency of the source

% Time scheme:
% 1= Explicit Newmark-alpha with alpha=1/2, beta=1/2, gamma=1
% 2= Centered difference = explicit Newmark with beta=0, gamma=1/2 (alpha=1)
if ~exist('TIME_SCHEME','var'), TIME_SCHEME=2; end

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
% so it's simple to cut windows of reflected phases for comparison
dt = 2*L/ceil(2*L/dt0);
disp(sprintf('Changed CFL from %f to %f',CFL,CFL*dt/dt0))
CFL = CFL*dt/dt0;
NT = ceil(TMAX/dt); % number of timesteps

half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;
if TIME_SCHEME==1
  alpha=1/2;
else
  alpha=1;
end

if ABSO_TOP
  BcTopC = sqrt(rho*mu);
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho*mu);
  M(1) = M(1)+half_dt*BcBottomC;
end

% treat implicitly the diagonal part of the damping matrix
%if NEL_VISC
%  Kdiag = zeros(NGLL,NEL_VISC);
%  for e=1:NEL_VISC,
%    ix = iglob(:,e);
%    Kdiag(:,e) = diag(K(:,:,e));
%    M(ix) = M(ix) + 0.1* 0.5*tvisc*Kdiag(:,e);
%  end
%end

d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);
f = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
if ~F_IS_WAVE, [Fdist,Fix] = min( abs(coor-Fx) ); end
% source time function 
Ft0 = 1.5/Ff0; % delay
t = (0:NT-1)'*dt+alpha*dt;
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
OUTa = zeros(OUTnx,NT);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Two time-integration schemes:
% 1= Explicit Newmark-alpha with alpha=1/2, beta=1/2, gamma=1
% 2= Centered difference = explicit Newmark with beta=0, gamma=1/2 (alpha=1)

for it=1:NT,

  if TIME_SCHEME==1
   % prediction of mid-step displacement:
   % d_mid = d_old + 0.5*dt*v_old
    d = d + half_dt*v; 
  else
   % update
    d = d + dt*v + half_dt_sq*a; 
   % prediction
    v = v + half_dt*a;
  end
  
 % internal forces at mid-step -K*d(t+1/2) 
  f(:) = 0;
  for e=1:NEL,
    ix = iglob(:,e);
    f(ix) = f(ix) - K*d(ix) ;
  end 

 % viscous layer:
  for e=1:NEL_VISC
    ix = iglob(:,e);
    f(ix) = f(ix) - tvisc* K*v(ix); 
  end

 % add external forces
  if ~F_IS_WAVE, f(Fix) = f(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, f(nglob) = f(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
       %boundary term = -impedance*v_out + impedance*v_in
       %              = -impedance*(v-v_in) + impedance*v_in
       %              = -impedance*(v-2*v_in)
      f(1) = f(1) -BcBottomC*( v(1) -2*Ft(it) );
    else
      f(1) = f(1) -BcBottomC*v(1);
    end
  end

 % acceleration: a = (-K*d +F)/M
  a = f ./M ;
  
  if TIME_SCHEME==1
   % update
   % v_new = v_old + dt*a_new;
   % d_new = d_old + dt*v_old + 0.5*dt^2*a_new
   %       = d_mid + 0.5*dt*v_new
    v = v + dt*a;
    d = d + half_dt*v;
  else
   % correction
    v = v + half_dt*a;
  end


%------------------------------------------
% STEP 4: OUTPUT
  
  OUTd(:,it) = d(OUTix);
  OUTv(:,it) = v(OUTix);
  OUTa(:,it) = a(OUTix);

end % of time loop
