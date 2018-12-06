% analytic_sh2d scalar 2D wave equation seismograms computed analytically
%               by direct integration of the Green's function 
% 		  G(r,t) = H(t-r/c) / (2*pi*sqrt(t^2-r^2/c^2))
% 		Solution:
% 		  u(r,t) = 1/(2*pi) integral_0_to_(t-ta) [ f(t')/sqrt((t-t')^2-ta^2) ] dt'
% 		with ta=r/c the arrival time
%
% Example:
% [u,stats] = analytic_sh2d(ta,f,t,tol)
%
% INPUT		ta	arrival time
%		f	source time function handle f(t)
%		t	times
%		tol	absolute tolerance for quadl 
%			(make sure it is much smaller than the maximum signal)
% 
% OUTPUT	u	seismogram
%		stats	structure containing algorithmic information:
%			n	total number of function calls
%            		np	number of function calls per non-null timestep
%
function [u,stats] = analytic_sh2d(ta,f,t,tol)

nt= length(t);
u = zeros(nt,1);

if ~exist('tol','var'), tol=max(abs(f(t)))/sqrt(ta)/pi*1e-6 , end

% the integrand has a sqrt singularity at t'=t-ta
% by a change of variable (see Numerical Recipes eq. 4.4.5):
% G(r,t) = 1/(2*pi) integral_0_to_sqrt(t-ta) [ 2x f(t-ta-x^2)/sqrt((t-(t-ta-x^2))^2-ta^2) ] dx
%        = 1/pi integral_0_to_sqrt(t-ta) [ f(t-ta-x^2)/sqrt(x^2+2*ta) ] dx
n=0;
np=0;
for it= 1:nt,
  if t(it)>ta
    dT = t(it)-ta;
    [u(it),ni] = quadl(@(x) f(dT-x.^2)./sqrt(x.^2+2*ta) ,0, sqrt(dT), tol ); 
    n=n+ni;
    np=np+1;
  end
end
u = u/pi;

stats.n = n;
stats.np = n/np;
