% Dispersion analysis for the central difference time-scheme
% =Newmark with alpha=1, beta=0, gamma=1/2
% See Hughes (1987) page 498

wdt=logspace(-2,log10(2),100); % omega*dt

A1 = 1-0.5*wdt.^2;
A2 = 1;
lp = A1 + sqrt(A1.^2-A2); % lambda_plus
lm = A1 - sqrt(A1.^2-A2); % lambda_minus

% lambda = exp(s*dt)
spw = log(lp)./wdt; % s / omega
smw = log(lm)./wdt;

%plot(wdt,imag(spw),'b-o', wdt,imag(smw),'r-')
% to make sure there is no attenuation, real(s) should be zero:
%plot(wdt,real(spw),'b--', wdt,real(smw),'r--')

loglog(wdt,abs(imag(spw))-1, wdt,abs(imag(smw))-1 ,wdt,(3/8-1/3)*wdt.^2,'--')
loglog(wdt,abs(imag(spw))-1, wdt,abs(imag(smw))-1 ,wdt,1/(3*8)*wdt.^2,'--')
xlabel('\omega \times \Delta t')
ylabel('relative frequency error \Delta\omega/\omega')
