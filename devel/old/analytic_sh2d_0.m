% Green's function for the scalar 2D wave equation:
% G(x,y,t) = H(t-r/c) / (2*pi*sqrt(t^2-r^2/c^2))  with r=sqrt(x^2+y^2)
%
% solution:
% u(r,t) = 1/(2*pi) integral_0_to_(t-r/c) [ f(t')/sqrt((t-t')^2-r^2/c^2) ] dt'

function [u,t] = analytic_sh2d_0(tar,f,dt,tmax)

nt=ceil(tmax/dt);
u = zeros(nt,1);
t = [1:nt]*dt;
sub = 10;
fs= f( (0:sub*nt)*dt/sub );
dtsub=dt/sub;

% start computing right after arrival time (it>tar/dt)
ita = ceil(tar/dt);

for it=ita:nt,
  ti = it*dt;
  tjm = 0;
  us = 0;
  tup = ti-tar;
  jtup = floor(tup/dtsub);
  if jtup>0
   for jt = 1:jtup,
    tj = jt*dtsub;
    tmid = 0.5*(tjm+tj);
    dint = (int0(tj,ti,tar)-int0(tjm,ti,tar))*(fs(jt)+fs(jt+1))*0.5 ...
         + (int1(tj,ti,tar,tmid)-int1(tjm,ti,tar,tmid))*(fs(jt+1)-fs(jt))/dtsub ;
    us = us+dint;
    tjm=tj;
   end
  else
    jt=0;
  end
  if tjm<tup
    dint = (int0(tup,ti,tar)-int0(tjm,ti,tar))*(fs(jt+1)+f(tup))*0.5 ;
    us = us+dint;
  end
  u(it) = us;
end
u = u/(2*pi);

%-------
function s = int0(x,a,b)
dum = sqrt((a-x)^2-b^2);
s = -log(a-x + dum);
%-------
function s = int1(x,a,b,c)
dum = sqrt((a-x)^2-b^2);
s = dum - (a-c)*log(a-x + dum) ;
