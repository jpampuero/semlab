% testing implicit EG-alpha scheme
% in a single d.o.f. oscillator
%   M*u_tt = - K*u + f(t)

M=1;
K=1;

CFL = 10;
NT = 1000;
dt = CFL/sqrt(K/M);
t=[1:NT]*dt;

r = 1.;
alpha_m = (2*r-1)/(r+1);
alpha_f = r/(1+r);
beta = 1/4*(1-alpha_m+alpha_f)^2;
gamma = 0.5-alpha_m+alpha_f;

D = (1-alpha_m)*M + ((1-alpha_f)*beta*dt^2)*K ;

% source time function (at mid-steps)
Ff0 = 0.015; % fundamental frequency
Ft0 = 1.5/Ff0; % delay
Ft = ricker( t'-0.5*dt, Ff0,Ft0);

d=0;
v=0;
a=0;

V=zeros(NT,1);

for it=1:NT,

  dlhs = d;
  d = d +dt*v + (dt^2*(0.5-beta))*a;
  dlhs = (1-alpha_f)*d + alpha_f*dlhs;
  v = v +(dt*(1-gamma))*a;

  lhs = -K*dlhs + Ft(it) - alpha_m*M*a;

  a = lhs/D ;

 % update
  v = v + (gamma*dt)*a;
  d = d + (beta*dt^2)*a;

  V(it)=v;

end
