fmax = 20;	% maximum frequency to be resolved
Npw = 1; 	% minimum number of elements per wavelength

[RHO,cs,eta,h]=textread('Oscar_Roma.tab','%*s%f%*f%f%f%f',...
                'commentstyle','matlab');
Q = 1./(2*eta);

cs2=[cs' ;cs'];cs2=cs2(:);
h(end) = cs(end)/fmax;
h2=[h' ;h'];h2=h2(:);
dep=cumsum(h);
dep2=[[0;dep(1:end-1)]';dep'];  dep2=dep2(:);

subplot(131)
plot(cs2,-dep2)
ylabel('z (m)')
xlabel('c_s (m/s)')
title('SH velocity model')
axis([0 1500 -inf 0])

subplot(132)
hol = h2./(cs2/fmax);
nel2 = ceil(Npw*hol);
plot( nel2,-dep2, '--',hol,-dep2 , nel2./hol, -dep2, ':')
xlabel('H / \lambda_{min} = H f_{max}/c_s')
legend('Elements per layer','layer thickness / \lambda_{min}','Elements per \lambda_{min}',4)
title(sprintf('%s\n%s',...
     'For mesh generation with f_{max} = 20 Hz',...
     '(resolution/dispersion criterion)'))
axis([0 7 -inf 0])


subplot(133)
dt = h2./nel2./cs2;
dt = dt/max(dt);
plot( dt, -dep2)
xlabel('\Delta t / \Delta t_{max}')
title(sprintf('%s\n%s',...
     'For time step selection',...
     '(stability criterion)'))
axis([0 1.2 -inf 0])

nel = nel2(1:2:end);
