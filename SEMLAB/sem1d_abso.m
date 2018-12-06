% SEM1D_ABSO 1D shear wave equation with absorbing boundaries
%
%   rho*u_tt = mu*u_xx + F(t)*delta(x-xs) 
%
% with a time dependent point source F(t) located at x=xs
% or an incident wave field F(t) from the bottom,
% zero initial conditions
% and stress free or absorbing boundary conditions.
%
% This script is intended for tutorial purposes.
%
% Jean-Paul Ampuero	jampuero@princeton.edu
%

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=10;		% domain size
NEL = 10;	% number of elements
P = 6;		% polynomial degree
CFL = 0.85; 	% stability number (<=0.85)
NT = 1000; 	% number of timesteps
ABSO_BOTTOM =0;	% absorbing boundary at x=0 ? if not: stress free boundary
ABSO_TOP = 0;	% absorbing boundary at x=L ? if not: stress free boundary
F_IS_WAVE = 0; 	% an incident wave (from bottom) or a point force ?
Fx = 5;		% location of point source (required if F_IS_WAVE=0)
Ff0 = 0.25; 	% fundamental frequency of the source
OUTx  = [0:L/NEL:L]'; % output receivers at these locations

%------------------------------------------
% STEP 2: INITIALIZATION

X = [0:NEL]'*L/NEL;
NGLL = P+1;
[iglob,coor,nglob] = mesh1d(X,NGLL);

%-- Set the physical properties of the medium
% Can be heterogeneous inside the elements
% and/or discontinuous across elements 
rho = ones(NGLL,NEL);		% density
mu  = ones(NGLL,NEL);		% shear modulus

[M,K] = BuildMK_1d(coor,iglob);

dt = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
half_dt = 0.5*dt;
t = (1:NT)' *dt;

%-- Absorbing boundary condition
% 	Traction = - impedance * velocity
% adds a term C*v to the equation at the boundary node
if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL)); % impedance
 % The mass matrix needs to be modified at the boundary
 % for the implicit treatment of the term C*v.
 % See the solver below for details.
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
  M(1) = M(1)+half_dt*BcBottomC;
end

d = zeros(nglob,1);
v = zeros(nglob,1);
a = zeros(nglob,1);

%-- External force (SOURCE TERM), a Ricker wavelet
% or velocity amplitude of an incoming wave
if ~F_IS_WAVE, [Fdist,Fix] = min( abs(coor-Fx) ); end
% source time function (at mid-steps)
Ft = src_timef( t-0.5*dt,'ricker', Ff0);

%-- Output arrays
OUTnx = length(OUTx);
% relocate the receivers to the nearest GLL node
for i=1:OUTnx, [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) ); end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,NT);
OUTv = zeros(OUTnx,NT);
OUTa = zeros(OUTnx,NT);

%------------------------------------------
% STEP 3: SOLVER  
% 
% M*a = -K*d -C*v +F
%
% Explicit Newmark (HHT-alpha) scheme with
% alpha=1/2, beta=1/2, gamma=1
%
% M*a(t+1) = -K*d(t+1/2) -C*v(t+1/2) +F(t+1/2)
% d(t+1) = d(t) + dt*v(t) +0.5*dt^2 *a(t+1)
% v(t+1) = v(t) + dt*a(t+1)
%

for it=1:NT,

%   1. predictor at mid-step t+0.5, assuming a(t+1)=0
%	dpre = d(t) +0.5*dt*v(t)
  d = d + half_dt*v; 

%   2. solve for a(t+1) in M*a(t+1) = -K*dpre +F(t+0.5)
%
  a(:) = 0;
  for e=1:NEL,
    ix = iglob(:,e);
    a(ix) = a(ix) - K(:,:,e)*d(ix) ;
  end 

 % add point forces
  if ~F_IS_WAVE, a(Fix) = a(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE 
   % incident wave, from bottom
   % boundary term = -impedance*v_outgoing + impedance*v_incoming
      a(1) = a(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
    else
      a(1) = a(1) -BcBottomC*v(1);
    end
  end

  a = a ./M ;

%   3. corrector
%	v(t+1) = v(t) + dt*a(t+1)
%	d(t+1) = dpre + 0.5*dt*v(t+1)
  v = v + dt*a;
  d = d + half_dt*v;


%------------------------------------------
% STEP 4: OUTPUT
  
  OUTd(:,it) = d(OUTix);
  OUTv(:,it) = v(OUTix);
  OUTa(:,it) = a(OUTix);

end % of time loop

PlotSeisTrace(OUTx,t,OUTv);
