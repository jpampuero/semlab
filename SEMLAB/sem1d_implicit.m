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
% THIS VERSION: implicit generalized alpha scheme
%		(J. Chung and G.M. Hulbert, JAM 1993)
%		which allows for a user-controlled dissipation at high-frequency 
%		while mantaining low dissipation at low-frequency.
%		As it is an implicit scheme, large timesteps are allowed.
%		However, it introduces dispersion. 
%		In this version the inversion of the dynamic matrix 
%		is done by brute force.

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=20;		% domain size
NEL = 20;	% number of elements
P = 6;		% polynomial degree
CFL = 4; 	% stability number (0.85 in explicit, can be larger in implicit)
NT = 125; 	% number of timesteps
ABSO_BOTTOM =0;	% absorbing boundary at x=0 ?
ABSO_TOP = 0;	% absorbing boundary at x=L ?

r = 1.; 	% dissipation parameter of the generalized-alpha scheme
		% 1 = no dissipation
		% 0 = maximal high-frequency dissipation

%------------------------------------------
% STEP 2: INITIALIZATION

X = [0:NEL]'*L/NEL;
NGLL = P+1; % number of GLL nodes per element
[iglob,coor,nglob] = mesh1d(X,NGLL);	

rho = ones(NGLL,NEL);		% density
mu  = ones(NGLL,NEL);		% shear modulus
[M,Kloc] = BuildMK_1d(coor,iglob,rho,mu);

% Assemble the global stiffness matrix
% Not required if inverting by iterative method instead 
K = zeros(nglob,nglob);
for e=1:NEL,
  ig = iglob(:,e);
  K(ig,ig) = K(ig,ig) + Kloc(:,:,e);
end 

alpha_m = (2*r-1)/(r+1);
alpha_f = r/(1+r);
beta = 1/4*(1-alpha_m+alpha_f)^2;
gamma = 0.5-alpha_m+alpha_f;

dt = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
half_dt = 0.5*dt;

D = (1-alpha_m)*diag(M) + ((1-alpha_f)*beta*dt^2)*K ;
%invD = inv( (1-alpha_m)*diag(M) + ((1-alpha_f)*beta*dt^2)*K );

if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL));
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  M(1) = M(1)+half_dt*BcBottomC;
end

d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);
lhs = zeros(nglob,1);
dlhs = zeros(nglob,1);

% External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
F_IS_WAVE = 0;
if ~F_IS_WAVE,
  Fx = L/2;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 0.25; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = src_timef( (1:NT)'*dt-0.5*dt, 'ricker',Ff0,Ft0);

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe:L]';
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
%

for it=1:NT,

  dlhs = d;
  d = d +dt*v + (dt^2*(0.5-beta))*a;
  dlhs = (1-alpha_f)*d + alpha_f*dlhs;
  v = v +(dt*(1-gamma))*a;

 % internal forces at mid-step -K*d(t+1/2) 
  lhs(:) = 0; % store -K*d in a global array
  for e=1:NEL,
    ix = iglob(:,e);
    lhs(ix) = lhs(ix) - Kloc(:,:,e)*dlhs(ix) ;
  end 
  lhs = lhs - alpha_m*M.*a;

 % add external forces
  if ~F_IS_WAVE, lhs(Fix) = lhs(Fix) + Ft(it); end

% % absorbing boundary
%  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
%  if ABSO_BOTTOM, 
%    if F_IS_WAVE % incident wave, from bottom
%     %boundary term = -impedance*v_outgoing + impedance*v_incoming
%      a(1) = a(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
%    else
%      a(1) = a(1) -BcBottomC*v(1);
%    end
%  end

 % acceleration: a = (-K*d +F)/M
  %a = invD*lhs ;
  a = D\lhs ;

 % update
  v = v + (gamma*dt)*a;
  d = d + (beta*dt^2)*a;


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
PlotSeisTrace(OUTx,t,OUTv);
