% SEM1D	applies the Spectral Element Method
% to solve the 1D SH wave equation, 
% with stress free or absorbing boundary conditions,
% zero initial conditions
% and a time dependent force source or incident plane wave.
%
% THIS VERSION for Matt Haney: linear slip fault in the middle
%
% Jean-Paul Ampuero	ampuero@erdw.ethz.ch
%

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=20;		% domain size
NEL = 20;	% number of elements, must be even to place the fault at the middle
P = 6;		% polynomial degree
CFL = 0.8; 	% stability number, smaller than without linear-slip-fault (0.85)
NT = 500; 	% number of timesteps
ABSO_BOTTOM =1;	% absorbing boundary at x=0 ?
ABSO_TOP = 1;	% absorbing boundary at x=L ?
LSFeta = 1; 	% fault compliance


%------------------------------------------
% STEP 2: INITIALIZATION

% Build the "macro-mesh"
% The interval [0,L] is divided into NEL non overlapping elements
% The elements are defined by their "control nodes" X
X = [0:NEL]'*L/NEL;

% Build the SEM grid
NGLL = P+1; % number of GLL nodes per element
[iglob,coor,nglob] = mesh1d(X,NGLL);	

% add a fault at the middle: introduce a split-node 
ix = iglob(NGLL,NEL/2); % global index of the middle node
coor = [coor(1:ix); coor(ix); coor(ix+1:end)] ; % insert a split-node
iglob(:,NEL/2+1:end) = iglob(:,NEL/2+1:end) +1; % update the local-global index map
nglob = nglob+1;
LSFiglob1 = iglob(NGLL,NEL/2);
LSFiglob2 = iglob(1,NEL/2+1);

% Set the physical properties of the medium
% Can be heterogeneous inside the elements
% and/or discontinuous across elements 
rho = ones(NGLL,NEL);		% density
mu  = ones(NGLL,NEL);		% shear modulus

% Mass and stiffness matrices
[M,K] = BuildMK_1d(coor,iglob,rho,mu);

dt = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
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
F_IS_WAVE = 1;
if ~F_IS_WAVE,
  Fx = 5; % in case of a point force
  [Fdist,Fix] = min( abs(coor-Fx) ); % relocate to nearest GLL node
end
Ff0 = 0.4; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
% source time function (at mid-steps)
Ft = ricker( (1:NT)'*dt-0.5*dt, Ff0,Ft0);
%Ft = ricker( (0:NT)'*dt, Ff0,Ft0);

% output arrays
OUTx  = X; % receivers locations
OUTnx = length(OUTx);
% relocate to nearest GLL node
OUTix = zeros(OUTnx,1);
OUTdist = zeros(OUTnx,1);
for i=1:OUTnx,
  [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) );
end
OUTd = zeros(OUTnx,NT);
OUTv = zeros(OUTnx,NT);
OUTa = zeros(OUTnx,NT);

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

 % linear slip fault
  tau = (d(LSFiglob2)-d(LSFiglob1))/LSFeta;
  a(LSFiglob1) = a(LSFiglob1) +tau;
  a(LSFiglob2) = a(LSFiglob2) -tau;

 % absorbing boundaries
  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
     %boundary term = -impedance*v_out + impedance*v_in
     %              = -impedance*(v-v_in) + impedance*v_in
     %              = -impedance*(v-2*v_in)
      a(1) = a(1) -BcBottomC*( v(1) -2*Ft(it) );
      %a(1) = a(1) -BcBottomC*( v(1) -(Ft(it)+Ft(it+1)) );
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
  
  OUTd(:,it) = d(OUTix);
  OUTv(:,it) = v(OUTix);
  OUTa(:,it) = a(OUTix);

end % of time loop

t = (1:NT) *dt;
%Ft(1)=[];

figure(1)
PlotSeisTrace(OUTx,t,OUTv);


%---------------------------------------------
% STEP 5: COMPARE TO ANALYTICAL SOLUTION


