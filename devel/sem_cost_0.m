% Computational cost = total number of multiplications in the computation of internal forces
% needed to propagate a wave over NL wavelengths

PROB = {'2DPSV', '3D'};
PROB = {'3D'}; % 1D, 2DSH, 2DPSV, 3D
NL = [10 50 100];
epsilon = 10.^(-[1:4]);
NL =1;
epsilon = 10.^(-[2:5]);

load P_dx_wmax
P=P(2:end)';
wmax = wmax(2:end)';
%wmax = 7/12 *P.*(P+2); % approximation 1
%wmax = 0.64*P.^2 +0.57*P +1; % approximation 2

A = (factorial(P)./factorial(2*P)).^2 ./(2*P.*(2*P+1));
%A = (exp(1)/4./P).^(2*P) ./(4*P.*(2*P+1)); % with Stirling's approximation
C = 2;

nj=length(PROB);
ni=length(NL);

for jplot=1:nj,

figure(1)
switch PROB{jplot}
  case '1D'   , E = P+1; D=1;
  case '2DSH' , E = 4*(P+2); D=2;
  case '2DPSV', E = 8*(P+3); D=2;
  case '3D'   , E = 9*(2*P+11); D=3;
end
E = E.*(1+P).^D;
Wmax = wmax*sqrt(D);
Cost1 = (2*pi)^(1+D).* E .* Wmax/C ;

% loop on NL and epsilon
for iplot=1:ni,
 for k=1:length(epsilon),
  Cost(:,k) = Cost1 .* NL(iplot).^((1+D)*(1+1./(2*P))) .* (A/epsilon(k)).^((1+D)./(2*P)) ;
  [cm(k,iplot),im(k,iplot)]=min(Cost(:,k));
 end
  subplot(ni,nj,nj*(iplot-1)+jplot)
  semilogy(P,Cost, P(im(:,iplot)),cm(:,iplot),'o')
  if iplot==1, title(PROB{jplot});  end
  if jplot==1, ylabel(['L/\lambda = ' num2str(NL(iplot))]); end
  if jplot==nj, ylabel('Total computational cost'), end
end
xlabel('Polynomial order')
%axis([0 20 1e10 1e12])

  figure(2)
  subplot(1,nj,jplot)
  semilogx(epsilon,P(im),'-o')
  xlabel('Accuracy target \epsilon')
  if jplot==1
    ylabel('Optimal order p_{opt}'); 
    legend(['L/\lambda = ' num2str(NL(1))],num2str(NL(2)),num2str(NL(3)))
  end
  title(PROB{jplot})
  
end

figure(1)
subplot(ni,nj,1)
legend('\epsilon = 10%','1%','0.1%','0.01%')
