% Dispersion analysis for the explicit (predictor-corrector) versions
% of the Newmark and HHT-alpha time-schemes

eta=0; alpha = 1; beta = 1; gamma = 1; % Raul HHT-alpha without KV
eta=0; alpha = 1/2; beta = 1/2; gamma = 1;
eta = 0.1; alpha=1; beta=0; gamma=1/2; % centered difference with Kelvin-Voigt damping

NPASS = 2;

% With the notation d = displacement, v = velocity*dt, a = acceleration*dt^2
% and w = omega*dt :

w=logspace(-2,log10(2.2),100);
w=linspace(0.01,5.2,100);

% Predictor:
%   dp = d[n] + v[n] + (1/2-beta)*a[n]
%   vp = v[n] + (1-gamma)*a[n]
%
% xp = PRED*x[n]  where x = (d,v,a)
PRED = [[1 1 1/2-beta];...
        [0 1 1-gamma] ; ...
	[0 0 0]];

% Solver:
%   a[n+1] = - w^2*( alpha*dp +(1-alpha)*d[n] + eta*(alpha*vp +(1-alpha)*v[n]) )
% Corrector:
%   d[n+1] = dp + beta*a[n+1]
%   v[n+1] = vp + gamma*a[n+1]
%
% LHS*x[n+1] = (ARP+ARO)*xp + ARN*x[n]
%
% For next pass:
% LHS*x[n+1]_new = ARP*xp + ARO*x[n+1]_old + ARN*x[n]
%
LHS = [[1 0 -beta] ; ...
       [0 1 -gamma]; ...
       [0 0 1]];
ARP = [[1 0 0] ; ...
       [0 1 0]; ...
       [0 0 0] ];
for k=1:length(w),
  ARO = [[0 0 0] ; ...
         [0 0 0]; ...
         [-alpha*w(k)^2 -eta*alpha*w(k)^2 0]];
  ARN = [[0 0 0] ; ...
         [0 0 0]; ...
         [-(1-alpha)*w(k)^2 -eta*(1-alpha)*w(k)^2 0]];
  % rhs matrix in LHS*x[n+1] = RHS*x[n]
  RHS = (ARP+ARO)*PRED + ARN;
  A = inv(LHS)*RHS;
  for ipass=2:NPASS,
    A = inv(LHS)*(ARP*PRED+ARO*A + ARN);
  end
  lambda = eig(A);
  rho(k) = max(abs(lambda));
  cxroots = lambda(imag(lambda)~=0);
  if ~isempty(cxroots)
    phi(k) = max(angle(cxroots));
  else
    phi(k) =0;
  end
end

subplot(221)
plot(w,rho)
hold all
ylabel('Amplification \rho')
subplot(223)
plot(w,phi./w)
hold all
ylabel('Phase velocity c_{num}')
xlabel('\omega \Delta t')

subplot(222)
loglog(w,abs(1-rho))
hold all
ylabel('|1-\rho|')
subplot(224)
loglog(w,abs(1-phi./w), w,w.^2/24,'--')
hold all
legend('this scheme','centered difference','Location','SE')
ylabel('|c_{num}-c|/c')
xlabel('\omega \Delta t')
