% HGLL = hgll(Z,ZGLL)
% Computes the value of the NGLL Lagrange polynomial interpolants through
% the NGLL Gauss-Lobatto Legendre points ZGLL(1:NGLL) at the point Z
% ZGLL and Z are within [-1,1]
%
function HGLL = hgll(Z,ZGLL)

% Adapted from gll.f90 in SEM2DPACK

  Z = Z(:);
  ZGLL= ZGLL(:);

  NGLL = length(ZGLL);
  N = NGLL - 1;
  ALFAN = N*(N+1);
%  HGLL = - (1-Z.*Z).*pndleg(Z,N) ./( ALFAN*pnleg(ZGLL,N).*(Z-ZGLL) );

 % avoid division by zero:
  HGLL = - (1-Z.*Z).*pndleg(Z,N) ./( ALFAN*pnleg(ZGLL,N) );
  %if abs(Z - ZGLL(I)) < 1e-5, HGLL(I) = 1; end  % EPS=1e-5
  mask = abs(Z - ZGLL) < 1e-5;
  HGLL = (~mask).*HGLL ./(Z-ZGLL+mask) +mask;

%--------------
% Compute the derivative of the Nth order Legendre polynomial at Z.
% Based on the recursion formula for the Legendre polynomials.
function PNDLEG = pndleg(Z,N)

  P1   = ones(N+1,1);
  P2   = Z;
  P1D  = zeros(N+1,1);
  P2D  = ones(N+1,1);
  P3D  = ones(N+1,1);
  for FK = 1:N-1,
   P3  = ((2*FK+1)*Z.*P2 - FK*P1)/(FK+1);
   P3D = ((2*FK+1)*P2 + (2*FK+1)*Z.*P2D - FK*P1D) /(FK+1);
   P1  = P2;
   P2  = P3;
   P1D = P2D;
   P2D = P3D;
  end
  PNDLEG = P3D;

%--------------
% Compute the value of the Nth order Legendre polynomial at Z.
% Based on the recursion formula for the Legendre polynomials.
function PNLEG = pnleg(Z,N)

  P1   = ones(N+1,1);
  P2   = Z;
  P3   = P2;
  for FK = 1:N-1,
   P3  = ((2*FK+1)*Z.*P2 - FK*P1)/(FK+1);
   P1  = P2;
   P2  = P3;
  end
  PNLEG = P3;
