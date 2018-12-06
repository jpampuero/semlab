% dispersion error for P=8
% assumes prior run of sem1d_dispersion with Ps=8;

% timestep / critical_timestep
gamma = 1/3; 
% gamma = 1;

% dimension
DIM = 2;

figure(3)
clf
% there are some annoying spikes (yet unknown origin)
% so I select some of the wavenumbers:
k= [4:4:length(rk)] ;
rks=rk(k); % k*h = wavenumber * elements size
rls = rks/2/pi; % element size / wavelength
ws=w(k)';

subplot(311)
plot(rks/2/pi,ws./rks)
axis([0 26/2/pi 0. 2.5])
ylabel('Phase velocity / wave speed')
%xlabel('k \times h = wavenumber \times element size')
%xlabel('h / \lambda = element size / wavelength')
xlabel('h / \lambda = wavelengths per element')
    % "NGLL" label
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    text(xl(2)/2,0.7*yl(2),sprintf('P = %u',P),'HorizontalAlignment','center')
    % two scales in the same plot
    v = get(gca,'position');
    npwl = [10 6 5 4 3 2];
    box off
    axes('position',v,'Ytick',[],'Xlim',xl/P, ... %(2*pi*P),...
         'xtick',1./npwl,'xticklabel',npwl,...
         'tickdir','in','XAxisLocation','top','color','none')
    xlabel('P \lambda / h = nodes per wavelength')
    axes('position',v,'ytick',[],'xtick',[],'color','none','box','on')

% compare to central difference dispersion error = (omega*dt)^2 /24
% with omega*dt = 2*pi/lambda*CFL*dxmin = CFL *k*h *dxmin/h
rlt=[0.1 4];
rkt=rlt*2*pi;
%ert = 1/24 *(0.85* rkt * (coor(2)-coor(1)) ).^2 ;
% theoretical time-discretization error (centered difference)
ert = 1/24 *(2*gamma/sqrt(DIM)/wmax_num(end)*rkt).^2 ;

% theoretical space-discretization error from Mulder (1999)
Ap = 0.5*(factorial(P)/factorial(2*P))^2 /(2*P+1)/P;
erx = Ap *rks.^(2*P);

subplot(5,1,[3:5])
%loglog(rks,abs((ws-rks)./rks),'-o',rks,0.1*(rks/19).^16,'--',rks,5*eps./rks,'--')
%axis([0.5 30 1e-20 100])
%loglog(rls,abs((ws-rks)./rks),'-',rls,0.1*(rks/19).^(2*P),'--',rls,50*eps./rks,'--',...
loglog(rls,abs((ws-rks)./rks),'-',rls,erx,'--', ...
       rlt,ert,'k--', rls,100*eps./rks,'--' )
axis([0.05 5 1e-15 100])
set(gca,'Ytick',logspace(-18,0,7))
grid on
grid minor
ylabel('\DeltaT / T = relative dispersion error ')
xlabel('h / \lambda = wavelengths per element')
legend('\DeltaT/T', ...
       '\DeltaT/T = A(p) (kh)^{2P}', ...
       '\DeltaT/T = B (C kh/\Omega_{max})^2',...
       '\DeltaT/T \sim round-off',...
       'Location','northwest')
    % two scales in the same plot
    xl = get(gca,'xlim');
    v = get(gca,'position');
    npwl = [100 50 20 10 6 5 4 3 2];
    box off
    axes('position',v,'Ytick',[],'Xlim',xl/P, ... %(2*pi*P),...
         'xtick',1./npwl,'xticklabel',npwl,'xscale','log',...
         'tickdir','in','XAxisLocation','top','color','none')
    xlabel('P \lambda / h = nodes per wavelength')
    axes('position',v,'ytick',[],'xtick',[],'color','none','box','on')
