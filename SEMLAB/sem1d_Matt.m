% SEM1D	applies the Spectral Element Method
% to solve the 1D SH wave equation, 
% with stress free or absorbing boundary conditions,
% zero initial conditions
% and a time dependent force source or incident plane wave.
%
% THIS VERSION for Matt Haney: linear slip fault in the middle.
% The fault is passive: it's a stiff interface across which a displacement discontinuity (slip) is allowed
% and slip is proportional to shear stress, as in Haney et al (2007, http://gji.oxfordjournals.org/content/170/2/933).
%
% Jean-Paul Ampuero	ampuero@gps.caltech.edu
%

%------------------------------------------
% STEP 1: SET PARAMETERS 

L=20;		% domain size
NEL = 20;	% number of elements, must be even to place the fault at the middle
P = 8;		% polynomial degree
CFL = 0.82/8; 	% stability number, smaller than without linear-slip-fault (0.85)
NT = 1024*8; 	% number of timesteps
ABSO_BOTTOM =1;	% absorbing boundary at x=0 ?
ABSO_TOP = 1;	% absorbing boundary at x=L ?
FAULT = 0;	% add a linear-slip fault
LSFeta = 1; 	% fault compliance
F_IS_WAVE = 1;	% source is incident wavefield, from bottom
OUTx  = [0:L/NEL:L]; 	% receivers locations


%------------------------------------------
% STEP 2: INITIALIZATION

X = [0:NEL]'*L/NEL;
NGLL = P+1;
[iglob,coor,nglob] = mesh1d(X,NGLL);	

% add a fault at the middle: introduce a split-node 
if FAULT
  ix = iglob(NGLL,NEL/2); % global index of the middle node
  coor = [coor(1:ix); coor(ix); coor(ix+1:end)] ; % insert a split-node
  iglob(:,NEL/2+1:end) = iglob(:,NEL/2+1:end) +1; % update the local-global index map
  nglob = nglob+1;
  LSFiglob1 = iglob(NGLL,NEL/2);
  LSFiglob2 = iglob(1,NEL/2+1);
end

rho = ones(NGLL,NEL);		% density
mu  = ones(NGLL,NEL);		% shear modulus
[M,K] = BuildMK_1d(coor,iglob,rho,mu);

dt = SetTimestep_1d(CFL, coor,iglob, sqrt(mu./rho));
half_dt = 0.5*dt;
half_dt_sq = 0.5*dt^2;
t = (1:NT)' *dt;

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

% Incident field, a Ricker wavelet
Ff0 = 0.2; 	% fundamental frequency
Ft0 = 1.5/Ff0; 	% delay
Ft = src_timef( t,'ricker', Ff0,Ft0);% source time function 
Fix = 1;

% output arrays
OUTnx = length(OUTx);
% relocate to nearest GLL node
for i=1:OUTnx, [OUTdist(i),OUTix(i)] = min( abs(coor-OUTx(i)) ); end
OUTx = coor(OUTix);
OUTd = zeros(OUTnx,NT);
OUTv = zeros(OUTnx,NT);
OUTa = zeros(OUTnx,NT);

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark scheme with
% alpha=1, beta=0, gamma=1/2
%

for it=1:NT,

 % update
  d = d + dt*v + half_dt_sq*a; 
 % prediction
  v = v + half_dt*a;
  a(:) = 0; 

 % internal forces -K*d(t+1) 
 % stored in global array 'a'
  for e=1:NEL,
    ix = iglob(:,e);
    a(ix) = a(ix) - K(:,:,e)*d(ix) ;
  end 

 % add external forces
  if ~F_IS_WAVE, a(Fix) = a(Fix) + Ft(it); end

 % linear slip fault
  if FAULT
    tau = (d(LSFiglob2)-d(LSFiglob1))/LSFeta;
    a(LSFiglob1) = a(LSFiglob1) +tau;
    a(LSFiglob2) = a(LSFiglob2) -tau;
  end

 % absorbing boundaries
  if ABSO_TOP, a(nglob) = a(nglob) -BcTopC*v(nglob); end
  if ABSO_BOTTOM, 
    if F_IS_WAVE % incident wave, from bottom
     %boundary term = -impedance*v_out + impedance*v_in
     %              = -impedance*(v-v_in) + impedance*v_in
     %              = -impedance*(v-2*v_in)
      a(1) = a(1) -BcBottomC*( v(1) -2*Ft(it) );
    else
      a(1) = a(1) -BcBottomC*v(1);
    end
  end

 % acceleration: a = (-K*d +F)/M
  a = a ./M ;

 % correction
  v = v + half_dt*a;

%------------------------------------------
% STEP 4: OUTPUT
  
  OUTd(:,it) = d(OUTix);
  OUTv(:,it) = v(OUTix);
  OUTa(:,it) = a(OUTix);

end % of time loop

subplot(221)
PlotSeisTrace(OUTx,t,OUTv);


%---------------------------------------------
% STEP 5: COMPARE TO ANALYTICAL SOLUTION

nbot = (OUTnx+1)/2;
va = zeros(OUTnx,NT);
for k=1:OUTnx,
  va(k,:) = specshift( OUTv(1,:)' ,-OUTx(k)/dt );
end

if FAULT
  wf = [ [0:NT/2] [-NT/2+1:-1] ]' *2*pi/(NT*dt) ;
  Z1 = sqrt(rho(NGLL,NEL/2)*mu(NGLL,NEL/2));
  Z2 = sqrt(rho(1,NEL/2+1)*mu(1,NEL/2+1));
  Dco = Z1+Z2 - 1i*wf*LSFeta*Z1*Z2 ;
  Rco = (Z1-Z2 - 1i*wf*LSFeta*Z1*Z2)./Dco ; % reflection coefficient
  Tco = 2*Z1./Dco; % transmission coefficient
  Rco = conj(Rco);
  Tco = conj(Tco);
  for k=1:OUTnx,
    if k<=nbot,
      vr = real(ifft( Rco.*fft(va(OUTnx+1-k,:)') ));
      va(k,:) = va(k,:) + vr'; 
    else
      vt = real(ifft( Tco.*fft(va(k,:)') ));
      va(k,:) = vt ; 
    end
  end
end

res = sqrt( sum((OUTv-va).^2,2)./sum(va.^2,2) ) ;

subplot(222)
PlotSeisTrace(OUTx,t,OUTv-va);
subplot(223)
plot(X,res,'-o')
xlabel('Distance')
ylabel('Normalized RMS')
subplot(224)
plot(X,res./X,'-o')
xlabel('Distance')
ylabel('Normalized RMS / Distance')
