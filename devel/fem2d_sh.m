% FEM2D_SH applies the Finite Element Method
% on a regular grid of equilateral triangles,
% introduced by Andrews (1973), 
% to the SH wave equation.
% Homogeneous medium.
% Periodic boundaries.

LX=30;
NELX = 240;
NELY = 160;
NEL = 2*NELX*NELY;
dx = LX/NELX;
dy = sqrt(3)/2*dx;
LY = NELY*dy;

% Physical properties of the medium
rho = 1;
mu  = 1;
vs = sqrt(mu./rho); 

CFL = 0.5;
dt = dx/vs *CFL;
NT = 600;
OUTdt = 20; % for output

% The mesh looks like a vertical zig-zag strip
% triangles A point down
% triangles B point up
% first element, bottom left, is A
% last element, top right, is B (NELY is even)
nglob = (NELX+1)*(NELY+1);
ngsiz = [NELX+1 NELY+1];
elsiz = [NELX*2 NELY];
iglob = zeros(3,NEL);
ElementType = zeros(NEL,1);

x = repmat([-NELX/2:NELX/2]'*dx,1,NELY+1);
x(:,2:2:NELY+1) = x(:,2:2:NELY+1)-0.5*dx;
y = repmat([0:NELY]*dy,NELX+1,1);
x = x(:);
y = y(:);
%Local to global numbering
ex=repmat([1:NELX]',NELY/2,1)';
%odd ey: A B ... A B
ey=repmat([1:2:NELY-1],NELX,1);
ey = ey(:)';
 %odd ex: Type A [SE NE NW]
  e = sub2ind(elsiz,2*ex-1,ey);
  ElementType(e) = 1; 
  iglob(1,e) = sub2ind(ngsiz,ex,ey);
  iglob(2,e) = sub2ind(ngsiz,ex+1,ey+1);
  iglob(3,e) = sub2ind(ngsiz,ex,ey+1);
 %even ex: Type B [SW SE NW]
  e = e+1;
  ElementType(e) = 2; 
  iglob(1,e) = sub2ind(ngsiz,ex,ey);
  iglob(2,e) = sub2ind(ngsiz,ex+1,ey);
  iglob(3,e) = sub2ind(ngsiz,ex+1,ey+1);

%even ey: B A ... B A
ey=repmat([2:2:NELY],NELX,1);
ey = ey(:)';
 %odd ex: Type B [SW SE NW]
  e = sub2ind(elsiz,2*ex-1,ey);
  ElementType(e) = 2; 
  iglob(1,e) = sub2ind(ngsiz,ex,ey);
  iglob(2,e) = sub2ind(ngsiz,ex+1,ey);
  iglob(3,e) = sub2ind(ngsiz,ex,ey+1);
 %even ex: Type A [SE NE NW]
  e = e+1;
  ElementType(e) = 1; 
  iglob(1,e) = sub2ind(ngsiz,ex+1,ey);
  iglob(2,e) = sub2ind(ngsiz,ex+1,ey+1);
  iglob(3,e) = sub2ind(ngsiz,ex,ey+1);

% boundary nodes
igLeft   = sub2ind(ngsiz,repmat(1,1,NELY+1),[1:NELY+1]);
igRight  = sub2ind(ngsiz,repmat(NELX+1,1,NELY+1),[1:NELY+1]);
igBottom = sub2ind(ngsiz,[1:NELX+1],repmat(1,1,NELX+1));
igTop    = sub2ind(ngsiz,[1:NELX+1],repmat(NELY+1,1,NELX+1));

%Type A: [SE NE NW]
cex(:,1) = 0.5*[0 1 -1]'/dx;
cey(:,1) = 0.5*[-1 0.5 0.5]'/dy;
%Type B: [SW SE NW]
cex(:,2) = 0.5*[-1 1 0]'/dx;
cey(:,2) = 0.5*[-0.5 -0.5 1]'/dy;

%Type A: [b d f]
csx(:,1) = [0 -1 1]'/(2*dx); 
csy(:,1) = [2 -1 -1]'/(2*sqrt(3)*dx);
%Type B: [a c e]
csx(:,2) = [1 -1 0]'/(2*dx);
csy(:,2) = [1 1 -2]'/(2*sqrt(3)*dx);

% nucleation parameters
NucDt = 200; % duration (in timesteps)
NucNel = 8;
nucstf = ones(NT+1,1);
nucstf(1:NucDt+1) = (1-cos(pi*[0:NucDt]/NucDt))/2;
nucstf = diff(nucstf)/dt;
NucEx1 =  NELX-NucNel+1;
NucEx2 =  NELX+NucNel;
NucEl = sub2ind(elsiz,[NucEx1:NucEx2]',repmat(NELY/2,2*NucNel,1));

% initialize fields
v = zeros(nglob,1);
f = zeros(nglob,1);
sx = zeros(1,NEL);
sy = zeros(1,NEL);
eply = zeros(1,NEL);

% for output seismograms
ix = [round(NELX/4):4:round(3*NELX/4)]';
OUTiglob = sub2ind(ngsiz,ix,repmat(round(NELY/4),length(ix),1) );
OUTxseis = x(OUTiglob); 
OUTyseis = y(OUTiglob); 
OUTnseis = length(OUTiglob);
OUTv = zeros(OUTnseis,NT);
 
%--------------------------------------------------------------------
% SOLVER

mu2dt = 2*mu*dt;
vloc=zeros(3,NEL);
cex=cex(:,ElementType);
cey=cey(:,ElementType);
dex=zeros(1,NEL);
dey=zeros(1,NEL);
csx=csx(:,ElementType);
csy=csy(:,ElementType);

%v = exp(-x.^2/2^2); % test plane wave

time = (1:NT)'*dt;

for it = 1:NT,

  eply(NucEl) = nucstf(it); % source is a plastic strain
  
 %velocity in local storage
  vloc = v(iglob);

 %strain
  dex = dot(cex,vloc);
  dey = dot(cey,vloc);

 %stress
  sx = sx + mu2dt*dex;
  sy = sy + mu2dt*(dey-eply);
  %sy = sy + mu2dt*dey;

 %internal forces in local storage
  floc = csx.*repmat(sx,3,1) + csy.*repmat(sy,3,1);

 %assemble internal forces
  f(:) = 0;
  for e=1:NEL,
    ig = iglob(:,e);
    f(ig) = f(ig) + floc(:,e);
  end

 %enforce periodicity
  f(igLeft) = f(igLeft)+f(igRight);
  f(igRight) = f(igLeft);
  f(igBottom) = f(igBottom) + f(igTop);
  f(igTop) = f(igBottom);

%  f((nglob-1)/2)=f((nglob-1)/2)+nucstf(it); %test point source

 %velocity update
  v = v + dt/rho *f;

 % output
  OUTv(:,it) = v(OUTiglob);
  if mod(it,OUTdt) == 0

    figure(1) % seismograms
    PlotSeisTrace(OUTxseis,time,OUTv);
    title('FEM2D SH velocity seismograms')

    figure(2)
    Plot2dSnapshot(x,y,v,iglob',[-1e-2 1e-2]);
    title('FEM2D SH velocity snapshot')
    hold on
    plot(OUTxseis,OUTyseis,'^')
    hold off
    drawnow
  end

end


