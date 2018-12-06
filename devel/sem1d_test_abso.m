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

%------------------------------------------
% STEP 1: MESH GENERATION
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X
NEL = 10;
L=500;
X = [0:NEL]'*L/NEL;

%------------------------------------------
% STEP 2: INITIALIZATION

P = 8; % polynomial degree
NGLL = P+1; % number of GLL nodes per element
NT = 1000; % number of timesteps

ABSO_TOP = 1;
ABSO_BOTTOM =0;

% The Gauss-Lobatto-Legendre points and weights
% and derivatives of the Lagrange polynomials H_ij = h'_i(xgll(j))
% were pre-tabulated for the usual range of NGLL.
% The xgll are in [-1,1]
[xgll,wgll,H] = GetGLL(NGLL);

iglob = zeros(NGLL,NEL);	% local to global index mapping
rho = zeros(NGLL,NEL);		% density
mu = zeros(NGLL,NEL);		% shear modulus
nglob = NEL*(NGLL-1) + 1;	% number of global nodes
coor = zeros(nglob,1);		% coordinates of GLL nodes
M = zeros(nglob,1);		% global mass matrix, diagonal
K = zeros(NGLL,NGLL,NEL);	% local stiffness matrix
CFL = 0.85; 			% stability number
dt = Inf;  			% timestep (set later)

for e=1:NEL, % FOR EACH ELEMENT ...

 % The table I = iglob(i,e) maps the local numbering of the 
 % computational nodes (i,e) to their global numbering I.
 % 'iglob' is used to build global data 
 % from local data (contributions from each element)
  iglob(:,e) = (e-1)*(NGLL-1)+(1:NGLL)';

 % Coordinates of the computational (GLL) nodes
  dxe = X(e+1)-X(e);
  coor(iglob(:,e)) = 0.5*(X(e)+X(e+1)) + 0.5*dxe*xgll ;

 % Physical properties of the medium
 % can be heterogeneous inside the elements
 % and/or discontinuous across elements 
  rho(:,e) = 2500;
  mu(:,e) = 2500*2887^2;

 % For this simple mesh the jacobian of the global-local
 % coordinate map is a constant
  dx_dxi = 0.5*dxe;

 % The diagonal mass matrix is stored in a global array. 
 % It is assembled here, from its local contributions
 % Nodes at the boundary between two elements get 
 % contributions from both.
  M(iglob(:,e)) = M(iglob(:,e)) + wgll .*rho(:,e) *dx_dxi;

% The stiffness matrix K is not assembled at this point
% (it is sparse, block-diagonal)
% We only store its local contributions
  W = mu(:,e).*wgll/dx_dxi;
%  K(:,:,e) = H * diag(W)* H';
  K(:,:,e) = H * ( repmat(W,1,NGLL).* H');

% The timestep dt is set by the stability condition
%   dt = CFL*min(dx/vs)
  vs = sqrt(mu(:,e)./rho(:,e)); 
  vs = max( vs(1:NGLL-1), vs(2:NGLL) );
  dx = abs(diff( coor(iglob(:,e)) ));
  dt = min(dt, min(dx./vs));
end %... of element loop
dt = CFL*dt;
half_dt = 0.5*dt;

% absorbing boundaries
% The mass matrix needs to be modified at the boundary
% for the implicit treatment of the term C*v.
if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL));
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  M(1) = M(1)+half_dt*BcBottomC;
end

% Initialize kinematic fields, stored in global arrays
d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
F_IS_WAVE = 0;
if ~F_IS_WAVE,
  Fx = 0;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 30; % fundamental frequency
Ft0 = 0.05; %1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = ricker( (1:NT)'*dt-0.5*dt, Ff0,Ft0);

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe/2:500]';
OUTnx = length(OUTx);
OUTix = zeros(OUTnx,1);
% relocate to nearest GLL node
OUTdist = zeros(OUTnx,1);
for i=1:OUTnx,
  [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) );
end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,OUTnt);
OUTv = zeros(OUTnx,OUTnt);
OUTa = zeros(OUTnx,OUTnt);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark-alpha scheme with
% alpha=1/2, beta=1/2, gamma=1
%

for it=1:NT,

 % prediction of mid-step displacement:
 % d_mid = d_old + 0.5*dt*v_old
  d = d + half_dt*v; 

 % internal forces at mid-step -K*d(t+1/2) 
  a(:) = 0; % store -K*d in a global array
  for e=1:NEL,
    ix = iglob(:,e);
    a(ix) = a(ix) - K(:,:,e)*d(ix) ;
  end 

 % add external forces
  if ~F_IS_WAVE, a(Fix) = a(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
     %boundary term = -impedance*v_outgoing + impedance*v_incoming
      a(1) = a(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
    else
      a(1) = a(1) -BcBottomC*v(1);
    end
  end

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
  
  if mod(it,OUTdt) == 0
    OUTit = OUTit+1;
    OUTd(:,OUTit) = d(OUTix);
    OUTv(:,OUTit) = v(OUTix);
    OUTa(:,OUTit) = a(OUTix);
  end

end % of time loop

figure(1)
t = (1:OUTit) *OUTdt*dt;
PlotSeisTrace(OUTx,t,OUTv);

TW = Ft0; 
V=2887;
p1 = floor( (Ft0 + OUTx/V - TW)/dt )+1;
p2 = floor( (Ft0 - (OUTx-2*L)/V -TW)/dt )+1;
tw = floor(2*TW/dt);

ista = 1;
u = OUTv(ista,:)';
nfft = 2^(nextpow2(tw)+1);
ferr = abs( fft(u(p2(ista):p2(ista)+tw,ista),nfft) ) ...
    ./ abs( fft(u(p1(ista):p1(ista)+tw,ista),nfft) );
ferr = ferr(1:nfft/2+1);
f = [0:nfft/2] /(nfft*dt);
k = find( f>Ff0*0.1 & f<Ff0*3 );
figure(2)
clf
%loglog(f,ferr, f(k),ferr(k),'k-o')
loglog( f(k),ferr(k),'k-o', ...
        f(k),2.5e-3*(f(k)/Ff0).^2, [Ff0 Ff0],[ferr(k(1)) ferr(k(end))],':')
grid on
ylabel('Relative error')
xlabel('Frequency (Hz)')
legend('err(f)','f^2',2)
err = norm(u(p2(ista):p2(ista)+tw,ista)) / norm(u(p1(ista):p1(ista)+tw,ista));
title( sprintf('Error = %f',err) )
