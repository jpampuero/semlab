% SRC_TIMEF Source time functions
%
% f = src_timef(t,name,arguments)
% f = src_timef(t,'ricker',f0 [,t0])
% f = src_timef(t,'gaussian',f0 [,t0])
% f = src_timef(t,'gabor',f0,t0,gamma [,theta])
% f = src_timef(t,'crack',tw,t0,tap_t0,tap_tw)
%
% INPUT		t	times
%		name	of a source time function:
%			'gaussian'
%			'ricker' second derivative of a Gaussian
%			'gabor'
%		f0	dominant frequency
%		t0	time offset [1.5/f0]
% 		gamma 	controls the width of the signal
%			and the spectral bandwidth
% 
% OUTPUT	f	source time function evaluated at times t
%

function f = src_timef(t,name,varargin)

fun = str2func(name);
f = fun(t,varargin{:});


%-- RICKER -------------------------------------------
function f = ricker(t,f0,t0)

% safe offset to guarantee almost zero initial value:
if nargin<3, t0=1.5/f0; end 

arg = pi*f0*(t-t0);
arg = arg.*arg;
f = (2*arg-1).*exp(-arg);

%-- GABOR wavelet -------------------------------------------
function f = gabor(t,f0,t0,gamma,theta)

if nargin<2, theta=pi/2; end 
if isempty(t0), t0 = 0.45*gamma/f0; end

w0 = 2*pi*f0;
f = exp(-( w0/gamma*(t-t0) ).^2) .* cos( w0*(t-t0)+theta);

%-- GAUSSIAN -------------------------------------------
function f = gaussian(t,f0,t0)

if nargin<3, t0=0.45/f0; end 

f = exp(-( 2*pi*f0*(t-t0) ).^2) ;

%-- CRACK smoothed and tapered -------------------------------
function f = crack(t,tw,t0, tap_t0,tap_tw)

if nargin<3, t0=0; end 
if isempty(tw), t0=0; end

t=t-t0;
f = 2/tw *(t>0).*( sqrt(t) - (t>tw).*sqrt(t-tw) );

if nargin>3,
  f = f.* exp( -((t>tap_t0).*(t-tap_t0)/tap_tw).^2 );
end
