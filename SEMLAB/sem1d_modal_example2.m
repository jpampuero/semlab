% We compute the maximum frequency (wmax) of the 1D wave equation
% discretized on a single spectral element.
% This is a useful information in setting stability criterion 
% for explicit time integration schemes:
%   wmax*critical_timestep = critical_CFL (= 2 for centered difference)

P  = (1:19);	% range of polynomial orders to analyze

for p = P,
  if p>1
    [eigvec,eigval,x] = sem1d_modal_analysis([0;1],p+1);
    wmax(p) = max(eigval);	% = maximum eigenvalue
    dx(p) = x(2);			% = smallest GLL spacing
  else
    wmax(p) = 2;
    dx(p) = 1;
  end
end
dx=dx(P); wmax=wmax(P);
disp(' P    dx_min     w_max')
disp(sprintf('%2u   %f   %f \n',[P; dx; wmax]));

figure(1)
loglog(dx,wmax,'o',dx,(7/3)./dx, '--')
xlabel('\Deltax_{min} / h')
ylabel('\omega_{max} \times h/c')
title('\omega_{max} \approx 7/3/\Deltax_{min}')

figure(2)
pw=polyfit(P,wmax,2);
loglog(P,wmax,'-o', P,polyval(pw,P))
xlabel('P')
ylabel('\omega_{max} \times h/c')
legend('numeric','quadratic fit')
