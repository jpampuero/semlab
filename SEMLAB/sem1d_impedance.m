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
% => analyze numerical impedance

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=10;		% domain size
NEL = 10;	% number of elements
P = 8;		% polynomial degree
CFL = 0.5; 	% stability number (0.85)
NT = 4096; 	% number of timesteps
ABSO_BOTTOM =0;	% absorbing boundary at x=0 ?
ABSO_TOP = 1;	% absorbing boundary at x=L ?

% add a viscous layer (one element) next to the boundary :
NEL_VISC = 0;
tvisc = 0.02; % viscous time /dt

%------------------------------------------
% STEP 2: INITIALIZATION

X = [0:NEL]'*L/NEL;
NGLL = P+1; % number of GLL nodes per element
[iglob,coor,nglob] = mesh1d(X,NGLL);	

rho = ones(NGLL,NEL);		% density
mu  = ones(NGLL,NEL);		% shear modulus
[M,K] = BuildMK_1d(coor,iglob,rho,mu);

dt = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
half_dt = 0.5*dt;

if ABSO_TOP
  BcTopC = sqrt(rho(NGLL,NEL)*mu(NGLL,NEL));
  M(nglob) = M(nglob)+half_dt*BcTopC;
end
if ABSO_BOTTOM
  BcBottomC = sqrt(rho(1,1)*mu(1,1));
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
F_IS_WAVE = 0;
if ~F_IS_WAVE,
  Fx = 0;
  [Fdist,Fix] = min( abs(coor-Fx) );
end
Ff0 = 3; % fundamental frequency
Ft0 = 2/Ff0; % delay
% source time function (at mid-steps)
Ft = src_timef( (1:NT)'*dt-0.5*dt,'gaussian', Ff0, Ft0);

% output arrays
OUTdt = 1;	% output every OUTdt timesteps
OUTit = 0;	% a counter
OUTnt = NT/OUTdt; % number of output timesteps
% output receivers at these locations
OUTx  = [0:dxe/2:L]';
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

  f(:) = 0;

 % internal forces at mid-step -K*d(t+1/2) 
  for e=1:NEL,
    ix = iglob(:,e);
    f(ix) = f(ix) - K(:,:,e)*d(ix) ;
  end 

 % viscous layer:
  for e=1:NEL_VISC
    ix = iglob(:,e);
    f(ix) = f(ix) - tvisc* K(:,:,e)*v(ix); 
  end

 % add external forces
  if ~F_IS_WAVE, f(Fix) = f(Fix) + Ft(it); end

 % absorbing boundary
  if ABSO_TOP, f(nglob) = f(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
     %boundary term = -impedance*v_outgoing + impedance*v_incoming
      f(1) = f(1) -BcBottomC*( v(1) -Ft(it) ) + BcBottomC*Ft(it);
    else
      f(1) = f(1) -BcBottomC*v(1);
    end
  end

 % acceleration: a = (-K*d +F)/M
  a = f ./M ;

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

t = (1:OUTit) *OUTdt*dt;
figure(1)
PlotSeisTrace(OUTx,t,OUTv);

figure(2)
ista=1;
subplot(311)
[freq,fv]=plot_spec(OUTv(ista,:),dt);
ylabel('velocity bottom boundary')
subplot(312)
[freq,fFt] = plot_spec(Ft,dt);
ylabel('source time function')
hold on; loglog( [freq(1) freq(end)], [1e-3 1e-3]*max(fFt),':');hold off
subplot(313)
loglog(freq,fFt./fv)
ylabel('impedance = f/v')
