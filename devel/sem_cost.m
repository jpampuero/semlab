% Computational cost = total number of multiplications in the computation of internal forces
% needed to propagate a wave over NL wavelengths

REPLOT = 0;
PROB = {'2DSH' '2DPSV', '3D'};
PROB = {'3D'}; % 1D, 2DSH, 2DPSV, 3D
epsilon = 10.^(-[3:5]); % = Delta_T / T
TIME_SCHEME = 'PEFRL'; 

load P_dx_wmax
P=P(2:end)';
wmax = wmax(2:end)';
%wmax = 7/12 *P.*(P+2); % approximation 1
%wmax = 0.64*P.^2 +0.57*P +1; % approximation 2

A = (factorial(P)./factorial(2*P)).^2 ./(2*P.*(2*P+1));
%A = (exp(1)/4./P).^(2*P) ./(4*P.*(2*P+1)); % with Stirling's approximation

% from time discretization:
switch TIME_SCHEME
  case 'IDEAL',	C=0; 	B=0; 	Q=0; 	N=0;
  case 'PV', 	C = 2; 	B = 1/24;Q = 2;N = 1;
  case 'CPV',  	C=2*sqrt(3); B=1/720; Q = 4; N=2;
  case 'PFR', 	C=1.5734; ; B = 6.61431e-2; Q = 4; N = 3;
  case 'PEFRL',	C=2.97633; B=1/12500; Q = 4; N=4;
end

nj=length(PROB);

if REPLOT
  for k=1:nj,figure(k),hold on,end
  linetype = '--';
else
  for k=1:nj, figure(k),clf, end
  linetype = '-';
end

for jplot=1:nj,

  figure(jplot)
  switch PROB{jplot}
    case '1D'   , E = P+1; D=1;
    case '2DSH' , E = 4*(P+2); D=2;
    case '2DPSV', E = 8*(P+3); D=2;
    case '3D'   , E = 9*(2*P+11); D=3;
  end
  E = E.*(1+P).^D;
  Wmax = wmax*sqrt(D);
  
  % loop on epsilon
  for k=1:length(epsilon),

   % optimal cost w.r.t. dt and h
    kh = (Q*D./(Q*D+2*P) *epsilon(k)./A ).^(1./(2*P)); 	% wavenumber*h
    dtw = (2*P./(Q*D+2*P) *epsilon(k)/B ).^(1/Q) ;	% omega*Delta_t
    Cost(:,k) = N*E ./( kh.^D .* dtw );
    ppw(:,k) = 2*pi*P./kh; 				% points per wavelength
    epw(:,k) = 2*pi./kh; 				% elements per wavelength
    dt(:,k) = dtw .* Wmax./kh/C ; 			% Delta_t / critical_Delta_t

   % optimal w.r.t. p
    [cm(k),im]=min(Cost(:,k));
    Pm(k) = P(im);
    ppwm(k) = ppw( im, k);
    epwm(k) = epw( im, k);
    dtm(k) = dt( im, k);

  end

  subplot(1,2,1)
  semilogy(P,Cost,linetype, Pm,cm,'ko') 
  axis([0 20 1e3 1e6])
  title(PROB{jplot})
  ylabel('Computational cost \Gamma')
  legend(num2str(epsilon','%g'))
  xlabel('Polynomial order')
  
  subplot(3,2,2)
  plot( P, ppw,linetype, Pm, ppwm,'ko' )
  axis([0 20 3 10])
  ylabel('Nodes per \lambda')

  subplot(3,2,4)
  plot( P, epw,linetype, Pm, epwm,'ko' )
  axis([0 20 0 5])
  ylabel('Elements per \lambda')

  subplot(3,2,6)
  plot( P, dt ,linetype, Pm,dtm,'ko')
  axis([0 20 0 1.1])
  ylabel('\Deltat / \Deltat_c')

  xlabel('Polynomial order')

end

%  figure(2)
%  semilogx(epsilon,Pm,[linetype 'o'])
%  xlabel('Accuracy \Delta T / T')
%  ylabel('Optimal order p_{opt}'); 
%  legend(PROB)

figure(5)
eta = 0.5; % Delta_t / Delta_t_critical
for G=[4.5 5 6 10],
  ex = A.*(2*pi*P/G).^(2*P);
  et = B*(eta*C./Wmax*2*pi.*P/G).^Q;
%  subplot(121)
  loglog(ex,et, '-s')
  hold all 
%  subplot(122)
%  loglog(ex+et,et./ex, '-s')
%  hold all 
end
%subplot(121)
e=[1e-6 1]; loglog(e,e,'k--',e,0.1*e,'k:',e,10*e,'k:', e,100*e,'k:')
hold off
axis([1e-6 1e-1 1e-4 1e-2])
xlabel('\epsilon_X'), ylabel('\epsilon_T')
legend('G = 4.5','5','6','10','\epsilon_T = \epsilon_X', 'location','NW')
text(0.02,6e-3,'p=2')
text(3e-5,1.2e-4,'p=19')


figure(6)
e = logspace(-6,-2,20);
ps = [2 3 4 6 10 15];
for p = ps,
  semilogx(e, 2*pi*p .*(A(find(P==p))./e).^(1/2/p) )
  hold all
end
hold off
axis([1e-6 1e-3 0 30])
legend(num2str(ps'))
xlabel('\epsilon_X')
ylabel('G')

figure(7)
ps = [2 4 6 10 15];
t = [0:2*pi/200:pi/2];
for p = ps,
  k = ( cos(t).^(2*p+2) + sin(t).^(2*p+2) ).^(-1/2/p) ;
  plot( cos(t).*k, sin(t).*k )
  hold all
end
hold off
axis equal square
legend(num2str(ps'), 'location','SW')
xlabel('\kappa_x / \kappa')
ylabel('\kappa_y / \kappa')
