% [alpha,beta] = GetRayleighPars(Q,F0,T0,NT,DT,WEIGHT,FIG)
%
% Computes the Rayleigh damping parameters alpha and beta
% (C = alpha*M + beta*K)
% aiming for a constant quality factor Q (1/2*eta) 
% within the Ricker source bandwidth,
% by minimizing
%   ( 1/Q - alpha/w - beta*w )^2
% for the discrete set of frequencies that are used
% in the spectral analysis (usually FFT) of the synthetic seismograms 
% and for which the spectrum of the Ricker is larger than 5% of its maximum.
% Optionnaly, the minimization can be weighted by the Ricker spectrum.
%
% INPUTS:	Q	target quality factor
%		F0,T0	central frequency and initial offset of the Ricker wavelet
%		NT,DT	number and size of timestep
%		WEIGHT	[0] flag for weighted least squares minimization
%		FIG	[0] flag for plotting results
%
function [a,b] = GetRayleighPars(Q,F0,T0,NT,DT,WEIGHT,FIG)

RickerSpectralThreshold = 0.05;
if ~exist('WEIGHT','var'), WEIGHT=0; end
if ~exist('FIG','var'), FIG=0; end

src = src_timef( (1:NT)'*DT,'ricker', F0,T0); 
ntfft = 2^nextpow2(NT);
fsrc = abs( fft(src,ntfft) );
fsrc = fsrc(1:ntfft/2+1);
k = find( fsrc >= RickerSpectralThreshold*max(fsrc) );
fsrcup = fsrc(k);
f = (0:ntfft/2)'/(ntfft*DT); 
omega = 2*pi*f(k);

if WEIGHT
  ab = [ fsrcup./omega fsrcup.*omega ] \ fsrcup/Q  ;
else
  ab = [ 1./omega omega ] \ ones(length(k),1)/Q  ;
end
a = ab(1);
b = ab(2);

if FIG
  f = f(2:end);
  fsrc = fsrc(2:end);
  loglog( f, fsrc/max(fsrc), f, a./(2*pi*f) + b*2*pi*f, [f(1) f(end)], [1 1]/Q )
  legend('Normalized source spectrum','Computed Q^{-1}','Target Q^{-1}',3)
  xlabel('Frequency')
end
