% reference seismogram
dt = 0.01;
Ff0=0.5;
Ft0=1.5/Ff0;
nt=ceil(2*Ft0/dt);
nt = 2^nextpow2(nt);
t= (0:nt-1)'*dt;
uref = src_timef(t,'ricker', Ff0,Ft0);

% seismogram with dispersion
omega = 2*pi* [0:nt/2 (-nt/2+1):-1]' /(nt*dt) ;
dispersion = 1/24 *(omega*dt).^2; % as in centered differences
travel_time = 100/Ff0;
u = fft(uref,nt) .* exp(1i*omega.*dispersion*travel_time);
u = real(ifft(u,nt));

plot(t,uref, t,u,'r')

% cross-correlation time delay estimates
npad=10; % zero-pad, should be > max expected time error
delay = xcorrshift(uref,[zeros(npad,1);u;zeros(npad,1)]);
delay = (delay +1-(npad+1)) *dt;

% quasi-analytical (Marquering et al 1999)
%vref = diff([0;uref])/dt;
%delay2 = sum( vref .*(uref-u) )/sum(vref.*vref)
vref = 1i*omega.*fft(uref) ;
delay2 = real( sum( vref.*conj(fft(uref-u)) )/sum(vref.*conj(vref)) );

omega_0 = norm(diff(diff([0;0;uref]))) /norm( diff([0;uref]) ) /dt;
delay3 = -1/24*(omega_0*dt)^2 *travel_time;

[delay delay2 delay3]
[delay delay2 delay3]/delay

