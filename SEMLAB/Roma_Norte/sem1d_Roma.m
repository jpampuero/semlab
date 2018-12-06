% This version: Roma Norte, Mexico D.F. (Oscar Ishizawa)
% Centered difference time stepping
% Rayleigh damping

% SEM1D	applies the Spectral Element Method
% to solve the 1D SH wave equation, 
% with Rayleigh damping,
% absorbing or free boundary conditions,
% zero initial conditions
% and a time dependent force source or incident plane wave.
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

P = 4; % polynomial degree
NT = 6000; % number of timesteps
CFL = 0.72; % stability number [0.85] 
            % NOTE: smaller CFL must be used for low Q (attenuation)
ABSO_TOP = 0;
ABSO_BOTTOM =1;
F_IS_WAVE = 1; % = 0 for point force, =1 for incident plane wave
Ff0 = 5; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
OUTdt = 20;	% output every OUTdt timesteps


%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X

%chose fmax and nb.elements in basement in Oscar_Roma.m
figure(1)
Oscar_Roma
NLAYERS = length(dep); %the model has 14 layers + basement
NEL = sum(nel);
% mesh from bottom to top
% and assign material numbers
MAT = []; X = [];
X(1) = -dep(end);
for k=NLAYERS:-1:1,
  X = [X; -dep(k)+(1:nel(k))'*h(k)/nel(k) ];
  MAT = [MAT; repmat(k,nel(k),1) ]; 
end

figure(2)
plot(h2./nel2,-dep2, X*0,X,'k-+')
axis([-2 60 -inf 0])
xlabel('c_s (m/s)')
ylabel('z (m)')

%------------------------------------------
% STEP 2: INITIALIZATION

NGLL = P+1; % number of GLL nodes per element

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);

iglob = zeros(NGLL,NEL);	% local to global index mapping
coor = zeros(NGLL,NEL);		% coordinates of GLL nodes
rho = zeros(NGLL,NEL);		% density
mu = zeros(NGLL,NEL);		% shear modulus
nglob = NEL*(NGLL-1) + 1;	% number of global nodes
M = zeros(nglob,1);		% global mass matrix, diagonal
Mbis = zeros(nglob,1);		% modified global mass matrix, diagonal
CM = zeros(nglob,1);		% M-proportional damping matrix, diagonal
K = zeros(NGLL,NGLL,NEL);	% local stiffness matrix
Mtmp = zeros(NGLL,NEL);		% local mass matrix, temporary

dt = Inf;  			% timestep (set later)

for e=1:NEL, % FOR EACH ELEMENT ...

 % The table I = iglob(i,e) maps the local numbering of the 
 % computational nodes (i,e) to their global numbering I.
 % 'iglob' is used to build global data 
 % from local data (contributions from each element)
  iglob(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';

 % Coordinates of the computational (GLL) nodes
  dxe = X(e+1)-X(e);
  coor(:,e) = 0.5*(X(e)+X(e+1)) + 0.5*dxe*xgll ;

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
  rho(:,e) = RHO(MAT(e));
  mu(:,e) = RHO(MAT(e))*cs(MAT(e))^2;

 % jacobian of the global-local coordinate map
  dx_dxi = 0.5*dxe;

 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
 % The local storage is kept temporarily to build the damping matrix
  Mtmp(:,e) = wgll .*rho(:,e) *dx_dxi;
  M(iglob(:,e)) = M(iglob(:,e)) + Mtmp(:,e);

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  W = mu(:,e).*wgll/dx_dxi;
  K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
  vs = sqrt(mu(:,e)./rho(:,e)); 
  vs = max( vs(1:NGLL-1), vs(2:NGLL) );
  dx = abs(diff( coor(:,e) ));
  dt = min(dt, min(dx./vs));
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;

% Rayleigh damping 
% parameters for Q=1, will be scaled later
figure(5)
[alphaC,betaC] = GetRayleighPars(1,Ff0,Ft0,NT,dt,0,1);
%alphaC=0;betaC=0;
% modified M for the implicit treatment of the diagonal part of the damping term
Mbis = M;
for e=1:NEL,
  ix = iglob(:,e);
  CM(ix) = CM(ix) + alphaC*Mtmp(:,e)/Q(MAT(e)) ;
  Mbis(ix) = Mbis(ix) +half_dt*betaC/Q(MAT(e))*diag(K(:,:,e)) ;
end
Mbis = Mbis + half_dt*CM;
clear Mtmp

% absorbing boundaries
% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL));
  Mbis(nglob) = Mbis(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  Mbis(1) = Mbis(1)+half_dt*BcBottomC;
end

% Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);


% External force (SOURCE TERM), a Ricker wavelet
if ~F_IS_WAVE,
  Fel=1; Fgll=3;	% located in element Fel, local node Fgll
  Fix = iglob(Fgll,Fel);
end
% source time function
Ft = src_timef( (1:NT)'*dt,'ricker', Ff0,Ft0);


% SEISMOGRAMS: 
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output these global nodes:
OUTnx = 2*NEL+1;
OUTix = zeros(OUTnx,1);
OUTx  = zeros(OUTnx,1);
kmid = round((NGLL+1)/2); % middle node
for e=1:NEL,
  OUTix(2*(e-1)+1) = iglob(1,e); 
  OUTx(2*(e-1)+1) = coor(1,e);
  OUTix(2*(e-1)+2) = iglob(kmid,e);
  OUTx(2*(e-1)+2) = coor(kmid,e);
end
OUTix(OUTnx) = iglob(NGLL,NEL);
OUTx(OUTnx) = coor(NGLL,NEL);
% output arrays
OUTd = zeros(OUTnx,OUTnt);
OUTv = zeros(OUTnx,OUTnt);
OUTa = zeros(OUTnx,OUTnt);

%------------------------------------------
% STEP 3: SOLVER  M*a + C*v + K*d = F
% Central difference
%

for it=1:NT,

 % update
  d = d + dt*v + half_dt_sq*a; 
 % prediction 
  v = v + half_dt*a;
  a(:) = 0; 

 % internal forces and K-proportional damping forces:
 %   -K*( d(t+1) +betaC*vpred(t+1) )
 % stored in global array 'a'
  for e=1:NEL,
    ix = iglob(:,e);
    a(ix) = a(ix) - K(:,:,e)*( d(ix) +betaC/Q(MAT(e))*v(ix) );
  end 

 % M-proportional damping forces
  a = a - CM.*v;  

 % add external forces
  if ~F_IS_WAVE, a(Fix) = a(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
      a(1) = a(1) -BcBottomC*( v(1) -Ft(it) );
    else
      a(1) = a(1) -BcBottomC*v(1);
    end
  end

 % acceleration: a = ( -C*v -K*d +F)/M
  a = a ./Mbis ;

 % correction
  v = v + half_dt*a;


%------------------------------------------
% STEP 4: OUTPUT
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;
    OUTd(:,OUTit) = d(OUTix);
    OUTv(:,OUTit) = v(OUTix);
    OUTa(:,OUTit) = a(OUTix);
  end

end % of time loop

t = (1:OUTit) *OUTdt*dt;
figure(3)
PlotSeisTrace(OUTx,t,OUTv);

figure(4)
subplot(211)
plot(t,OUTv(end,:),t,Ft(1:OUTdt:NT))
xlabel('Time (s)')
ylabel('Ground velocity (m/s)')
title('Surface seismogram')
legend('Seismogram','Incident wave',0)
%spectrum
subplot(212)
ntfft = 2^nextpow2(NT);
fsis = fft(OUTv(end,:)',ntfft);
fsis = abs( fsis(2:ntfft/2+1) );
fsrc = fft(Ft(1:OUTdt:NT),ntfft);
fsrc = 2*abs( fsrc(2:ntfft/2+1) ); % 2* for free surface
f = (1:ntfft/2)'/(ntfft*dt);
loglog(f,fsis,f,fsrc,f,fsis./fsrc)
xlabel('Frequency (Hz)')
ylabel('Spectral amplitude')
legend('Seismogram','Rock response','Amplification ratio',3)
