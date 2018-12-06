% f = gabor(t,f0,t0,gamma,theta)
% gamma controls the width of the signal, and the spectral bandwidth
function f = gabor(t,f0,t0,gamma,theta)
w0 = 2*pi*f0;
if isempty(t0), t0 = 0.45*gamma/f0; end
if ~exist('theta','var'), theta=[]; end
if isempty(theta), theta = pi/2; end
f = exp(-( w0/gamma*(t-t0) ).^2) .* cos( w0*(t-t0)+theta);
