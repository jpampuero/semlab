% At which value of G (nodes per wavelength) are the two components 
% of dispersion error (spatial and temporal) comparable ?
% How large is the dispersion error then ?

% timestep / critical_timestep
gammas = [1/2 1/3 1/5];

load P_dx_wmax

% select polynomial orders
P = [3:10];
Omega = wmax(P);

% time scheme parameters (accuracy and stability)
q=2; B = 1/24; C=2; % centered diff

% dimension
DIM = 2;

Ap = 0.5*(factorial(P)./factorial(2*P)).^2 ./(2*P+1)./P;

for g=gammas,
  G = 2*pi*P./ (B./Ap .*(g*C/sqrt(DIM)./Omega).^q).^(1./(2*P-q));
  err = Ap.* (2*pi*P./G).^(2*P);
  semilogy(G,err,'-s')
  hold all
end
hold off
legend(num2str(gammas'))
xlabel('G')
ylabel('\epsilon')
text(G(end),err(end)*1.5,num2str(P(end)))
text(G(1),err(1)*1.5,num2str(P(1)))


