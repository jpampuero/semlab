% An essential property of the Spectral Element Method is
% the "spectral convergence" which implies a very low dispersion error.
% The relative travel time misfit of a wave discretized by the SEM scales as 
% 	(lambda/h)^(-2*P) 
% where P is the polynomial order, lambda the wavelength and h the element size.
% The ratio lambda/h is the number of elements per wavelength.
% This is the expected convergence rate as a function of h
% for a P-th order method. The dispersion error is sometimes written as
% 	exp[ -2*P*log(lambda/h) ]
% to highlight an "exponential convergence rate" as a function of P.
% 
% To illustrate this spectral convergence property
% here we compute the eigenfrequencies of the wave equation
% on the unit 1D segment, discretized by the SEM
% with a range of polynomial orders P.
% The dispersion error is equivalent to the misfit of the discrete 
% eigenfrequencies with respect to their analytical values.

BC = 'D'; 	% boundary condition: Dirichlet or Neumann
NEL= 30;	% number of elements
P  = (2:9);	% range of polynomial orders to analyze

xmesh=(0:NEL)' /NEL;
figure(1)
clf

for p = P,

  % discrete frequencies from SEM:
  [eigvec,eigval,x,err]=sem1d_modal_analysis(xmesh,p+1,[BC BC]);
  if BC=='N', eigval=eigval(2:end); end % if Neumann skip zero frequency

  % analytical frequencies:
  k= pi*(1:length(eigval))' ;

  % relative frequency error (dispersion error):
  dval= abs(eigval-k)./k;

  % select frequencies with 
  % dispersion error > intrinsic error of the eigenvalue decomposition
  delta = err./k; % absolute error on eigenvalues
  isel = find(dval > 2*delta./k);

  npw = 2*NEL./isel; % = elements per wavelength
  loglog(npw,dval(isel)); hold on

end

hold off
title('Spectral convergence as (\lambda/h)^{-2p}')
xlabel('\lambda / h = number of elements per wavelength')
ylabel('Dispersion error \Deltaf / f')
grid on
